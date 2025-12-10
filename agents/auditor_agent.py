"""
Auditor Agent - 化学验证和审核智能体
负责验证分子合法性、毒性和类药性
"""

try:
    from agentscope.agents import AgentBase
    from agentscope.message import Msg
except ImportError:
    try:
        from agentscope import AgentBase, Msg
    except ImportError:
        from agents.orchestrator_agent import AgentBase, Msg
from typing import Optional, Union, List, Dict, Any
import json


class AuditorAgent(AgentBase):
    """审核智能体 - 统一验证化学合法性、毒性和类药性"""
    
    def __init__(
        self,
        name: str = "Auditor",
        **kwargs
    ):
        """初始化审核智能体
        
        Args:
            name: 智能体名称
            **kwargs: 其他参数
        """
        super().__init__(name=name, **kwargs)
        
        # 导入验证器
        self.validators = self._init_validators()
    
    def _init_validators(self) -> dict:
        """初始化各类验证器"""
        validators = {}
        
        try:
            from evaluators.chem_validator import ChemValidator
            validators['chem'] = ChemValidator()
        except Exception as e:
            print(f"警告: 化学验证器初始化失败: {e}")
            validators['chem'] = None
        
        try:
            from evaluators.similarity_checker import SimilarityChecker
            validators['similarity'] = SimilarityChecker()
        except Exception as e:
            print(f"警告: 相似性检查器初始化失败: {e}")
            validators['similarity'] = None
        
        return validators
    
    def reply(self, x: Optional[Union[Msg, List[Msg]]] = None) -> Msg:
        """审核各模型的输出结果
        
        Args:
            x: 输入消息，包含各模型的结果
        
        Returns:
            审核结果
        """
        if isinstance(x, list):
            x = x[-1]
        
        # 解析输入
        try:
            data = json.loads(x.content) if isinstance(x.content, str) else x.content
        except json.JSONDecodeError:
            return Msg(
                name=self.name,
                content=json.dumps({'success': False, 'error': 'Invalid input format'}),
                role="assistant"
            )
        
        model_results = data.get('model_results', {})
        task_config = data.get('task_config', {})
        
        # 审核每个模型的结果
        audit_results = {}
        
        for model_name, result in model_results.items():
            if not result.get('success'):
                audit_results[model_name] = {
                    'valid': False,
                    'error': result.get('error', 'Unknown error')
                }
                continue
            
            # 执行审核
            audit_result = self._audit_single_result(
                model_name,
                result,
                task_config
            )
            
            audit_results[model_name] = audit_result
        
        # 生成综合评估
        summary = self._generate_summary(audit_results, task_config)
        
        final_result = {
            'success': True,
            'audit_results': audit_results,
            'summary': summary
        }
        
        return Msg(
            name=self.name,
            content=json.dumps(final_result, ensure_ascii=False),
            role="assistant"
        )
    
    def _audit_single_result(
        self,
        model_name: str,
        result: Dict[str, Any],
        task_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """审核单个模型的结果
        
        Args:
            model_name: 模型名称
            result: 模型结果
            task_config: 任务配置
        
        Returns:
            审核结果
        """
        audit = {
            'model': model_name,
            'valid': True,
            'checks': {},
            'warnings': [],
            'errors': []
        }
        
        # 解析模型输出
        try:
            response_data = result.get('response', '')
            if isinstance(response_data, str):
                response_data = json.loads(response_data)
            
            output = response_data.get('output', '')
            if isinstance(output, str):
                # 尝试从文本中提取 JSON
                import re
                json_match = re.search(r'\{.*\}', output, re.DOTALL)
                if json_match:
                    output = json.loads(json_match.group())
                else:
                    output = {'raw': output}
        
        except Exception as e:
            audit['valid'] = False
            audit['errors'].append(f"解析输出失败: {str(e)}")
            return audit
        
        # 获取 SMILES (如果是分子生成任务)
        smiles = output.get('smiles', '') if isinstance(output, dict) else ''
        
        # 获取约束条件（提前获取，避免变量未定义错误）
        constraints = task_config.get('constraints', {})
        
        if smiles and self.validators.get('chem'):
            # 化学合法性验证
            chem_result = self.validators['chem'].validate_smiles(smiles)
            audit['checks']['chemistry'] = chem_result
            
            if not chem_result.get('valid'):
                audit['valid'] = False
                audit['errors'].append("分子结构不合法")
            
            # 类药性检查
            if chem_result.get('valid'):
                druglikeness = self.validators['chem'].check_druglikeness(smiles)
                audit['checks']['druglikeness'] = druglikeness
                
                # 分子量检查
                if 'max_mw' in constraints:
                    max_mw = constraints['max_mw']
                    actual_mw = druglikeness.get('properties', {}).get('molecular_weight', 0)
                    if actual_mw > max_mw:
                        audit['warnings'].append(
                            f"分子量 {actual_mw:.1f} 超过上限 {max_mw}"
                        )
                
                # QED 评分检查
                if 'min_qed' in constraints:
                    min_qed = constraints['min_qed']
                    actual_qed = druglikeness.get('qed_score', 0)
                    if actual_qed < min_qed:
                        audit['warnings'].append(
                            f"QED 评分 {actual_qed:.3f} 低于下限 {min_qed}"
                        )
            
            # 毒性检查
            if constraints.get('check_toxicity') and chem_result.get('valid'):
                toxicity = self.validators['chem'].check_toxicity(smiles)
                audit['checks']['toxicity'] = toxicity
                
                if toxicity.get('alerts'):
                    audit['warnings'].append(
                        f"发现 {len(toxicity['alerts'])} 个毒性警告"
                    )
        
        return audit
    
    def _generate_summary(
        self,
        audit_results: Dict[str, Any],
        task_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """生成审核摘要
        
        Args:
            audit_results: 各模型的审核结果
            task_config: 任务配置
        
        Returns:
            审核摘要
        """
        total_models = len(audit_results)
        valid_models = sum(1 for r in audit_results.values() if r.get('valid'))
        
        summary = {
            'total_models': total_models,
            'valid_models': valid_models,
            'pass_rate': valid_models / total_models if total_models > 0 else 0,
            'recommendations': []
        }
        
        # 生成建议
        if valid_models == 0:
            summary['recommendations'].append(
                "⚠️ 所有模型的输出都未通过验证，建议检查任务描述或放宽约束条件"
            )
        elif valid_models < total_models:
            summary['recommendations'].append(
                f"ℹ️ {total_models - valid_models} 个模型的输出未通过验证"
            )
        else:
            summary['recommendations'].append(
                "✓ 所有模型的输出都通过了基本验证"
            )
        
        # 统计警告
        total_warnings = sum(
            len(r.get('warnings', [])) for r in audit_results.values()
        )
        
        if total_warnings > 0:
            summary['recommendations'].append(
                f"ℹ️ 共发现 {total_warnings} 个警告，建议仔细审查"
            )
        
        return summary
