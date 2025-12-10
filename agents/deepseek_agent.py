"""
DeepSeek Agent - DeepSeek-V3 模型智能体
"""

try:
    from agentscope.agents import AgentBase
    from agentscope.message import Msg
except ImportError:
    try:
        from agentscope import AgentBase, Msg
    except ImportError:
        from agents.orchestrator_agent import AgentBase, Msg
from typing import Optional, Union, List
import json
from openai import OpenAI
from config import get_config


class DeepSeekAgent(AgentBase):
    """DeepSeek 模型智能体"""
    
    def __init__(
        self,
        name: str = "DeepSeek",
        **kwargs
    ):
        """初始化 DeepSeek 智能体
        
        Args:
            name: 智能体名称
            **kwargs: 其他参数
        """
        super().__init__(name=name, **kwargs)
        
        # 获取配置
        config = get_config()
        deepseek_config = config.get_deepseek_config()
        
        # 初始化 OpenAI 兼容客户端
        self.client = OpenAI(
            api_key=deepseek_config.get('api_key'),
            base_url=deepseek_config.get('base_url')
        )
        
        self.model = deepseek_config.get('model', 'deepseek-chat')
        
        # 加载 Prompt 模板
        self.prompts = self._load_prompts()
    
    def _load_prompts(self) -> dict:
        """加载 Prompt 模板"""
        from pathlib import Path
        import yaml
        
        prompts = {}
        prompts_dir = Path(__file__).parent.parent / "prompts"
        
        # 加载各类任务的 prompt
        prompt_files = {
            'molecule_generation': 'molecule_generation.yaml',
            'literature_extraction': 'literature_extraction.yaml',
            'admet_prediction': 'admet_protocol.yaml'
        }
        
        for task_type, filename in prompt_files.items():
            filepath = prompts_dir / filename
            if filepath.exists():
                with open(filepath, 'r', encoding='utf-8') as f:
                    prompts[task_type] = yaml.safe_load(f)
        
        return prompts
    
    def reply(self, x: Optional[Union[Msg, List[Msg]]] = None) -> Msg:
        """处理消息并返回结果
        
        Args:
            x: 输入消息
        
        Returns:
            处理结果
        """
        if isinstance(x, list):
            x = x[-1]
        
        # 解析任务配置
        try:
            task_config = json.loads(x.content) if isinstance(x.content, str) else x.content
        except json.JSONDecodeError:
            task_config = {
                'task_type': 'molecule_generation',
                'target': x.content
            }
        
        task_type = task_config.get('task_type', 'molecule_generation')
        target = task_config.get('target', '')
        
        # 构造 Prompt
        system_prompt = self._get_system_prompt(task_type)
        user_prompt = self._build_user_prompt(task_type, task_config)
        
        # 调用 DeepSeek API
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt}
                ],
                temperature=0.7,
                max_tokens=2000
            )
            
            content = response.choices[0].message.content
            
            # 封装结果
            result = {
                'model': 'DeepSeek-V3',
                'task_type': task_type,
                'target': target,
                'output': content,
                'metadata': {
                    'usage': {
                        'prompt_tokens': response.usage.prompt_tokens,
                        'completion_tokens': response.usage.completion_tokens,
                        'total_tokens': response.usage.total_tokens
                    }
                }
            }
            
            return Msg(
                name=self.name,
                content=json.dumps(result, ensure_ascii=False),
                role="assistant"
            )
            
        except Exception as e:
            error_result = {
                'model': 'DeepSeek-V3',
                'task_type': task_type,
                'error': str(e)
            }
            
            return Msg(
                name=self.name,
                content=json.dumps(error_result, ensure_ascii=False),
                role="assistant"
            )
    
    def _get_system_prompt(self, task_type: str) -> str:
        """获取系统 Prompt
        
        Args:
            task_type: 任务类型
        
        Returns:
            系统 Prompt
        """
        if task_type in self.prompts:
            return self.prompts[task_type].get('system_prompt', self._default_system_prompt())
        
        return self._default_system_prompt()
    
    def _default_system_prompt(self) -> str:
        """默认系统 Prompt"""
        return """你是一个专业的医药化学AI助手，精通分子设计、药物化学和ADMET预测。
你的任务是根据用户需求，提供准确、专业的药物化学建议和分析。

请遵循以下原则:
1. 所有分子结构使用SMILES格式表示
2. 提供详细的化学性质分析
3. 考虑类药性(drug-likeness)和成药性
4. 关注安全性和毒性问题
5. 提供可操作的建议

请以JSON格式返回结果。"""
    
    def _build_user_prompt(self, task_type: str, task_config: dict) -> str:
        """构建用户 Prompt
        
        Args:
            task_type: 任务类型
            task_config: 任务配置
        
        Returns:
            用户 Prompt
        """
        target = task_config.get('target', '')
        
        if task_type == 'molecule_generation':
            constraints = task_config.get('constraints', {})
            prompt = f"""请设计满足以下要求的分子:

目标: {target}

约束条件:
- 最大分子量: {constraints.get('max_mw', 500)}
- 最小QED评分: {constraints.get('min_qed', 0.5)}
- 检查毒性: {'是' if constraints.get('check_toxicity') else '否'}

请返回JSON格式，包含以下字段:
- smiles: 分子的SMILES表示
- name: 分子名称
- properties: 分子性质(MW, LogP, HBD, HBA, TPSA等)
- qed_score: QED评分
- rationale: 设计理由
"""
        
        elif task_type == 'literature_extraction':
            prompt = f"""请从医药文献中提取以下信息:

目标: {target}

请返回JSON格式，包含:
- key_findings: 关键发现
- molecules: 提及的分子列表(含SMILES)
- clinical_data: 临床数据
- references: 参考文献
"""
        
        elif task_type == 'admet_prediction':
            prompt = f"""请预测以下分子的ADMET性质:

目标: {target}

请返回JSON格式，包含:
- absorption: 吸收性质
- distribution: 分布性质
- metabolism: 代谢性质
- excretion: 排泄性质
- toxicity: 毒性预测
"""
        
        else:
            prompt = target
        
        return prompt
