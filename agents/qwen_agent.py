"""
Qwen Agent - 通义千问模型智能体
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


class QwenAgent(AgentBase):
    """Qwen 模型智能体 - 通过 DashScope 兼容接口调用"""
    
    def __init__(
        self,
        name: str = "Qwen",
        **kwargs
    ):
        """初始化 Qwen 智能体
        
        Args:
            name: 智能体名称
            **kwargs: 其他参数
        """
        super().__init__(name=name, **kwargs)
        
        # 获取配置
        config = get_config()
        qwen_config = config.get_qwen_config()
        
        # 初始化 OpenAI 兼容客户端 (DashScope 支持 OpenAI 格式)
        self.client = OpenAI(
            api_key=qwen_config.get('api_key'),
            base_url=qwen_config.get('base_url')
        )
        
        self.model = qwen_config.get('model', 'qwen-max')
        
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
        
        # 调用 Qwen API
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
                'model': 'Qwen-Max',
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
                'model': 'Qwen-Max',
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
        return """你是通义千问，一个专业的医药化学AI助手。你精通药物设计、化学信息学和生物医学。

你的专长包括:
1. 基于靶点的分子设计
2. ADMET性质预测
3. 药物文献分析
4. 化学结构优化
5. 成药性评估

请遵循以下原则:
- 使用SMILES格式表示分子结构
- 提供详细的科学依据
- 考虑实际可合成性
- 关注安全性和有效性
- 以JSON格式返回结构化数据

请专业、准确地完成任务。"""
    
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
            prompt = f"""分子设计任务:

【目标】{target}

【约束条件】
- 分子量上限: {constraints.get('max_mw', 500)} Da
- QED评分下限: {constraints.get('min_qed', 0.5)}
- 毒性检查: {'需要' if constraints.get('check_toxicity') else '不需要'}

【输出要求】
请以JSON格式返回，包含:
{{
  "smiles": "分子SMILES表示",
  "name": "建议的分子名称",
  "properties": {{
    "molecular_weight": 数值,
    "logp": 数值,
    "hbd": 数值,
    "hba": 数值,
    "tpsa": 数值
  }},
  "qed_score": 数值,
  "design_rationale": "设计思路说明",
  "synthesis_feasibility": "可合成性评估"
}}
"""
        
        elif task_type == 'literature_extraction':
            prompt = f"""文献信息提取任务:

【检索目标】{target}

【提取要求】
请以JSON格式返回，包含:
{{
  "key_findings": ["关键发现1", "关键发现2"],
  "molecules": [
    {{"name": "分子名", "smiles": "SMILES", "activity": "活性描述"}}
  ],
  "clinical_data": {{
    "phase": "临床阶段",
    "results": "结果摘要"
  }},
  "references": ["DOI或PMID"]
}}
"""
        
        elif task_type == 'admet_prediction':
            prompt = f"""ADMET性质预测任务:

【分子信息】{target}

【预测要求】
请以JSON格式返回ADMET分析:
{{
  "absorption": {{
    "caco2_permeability": "评估",
    "bioavailability": "预测值"
  }},
  "distribution": {{
    "vd": "分布体积",
    "bbb_penetration": "血脑屏障通透性"
  }},
  "metabolism": {{
    "cyp_substrates": ["CYP酶"],
    "metabolic_stability": "稳定性"
  }},
  "excretion": {{
    "clearance": "清除率",
    "half_life": "半衰期"
  }},
  "toxicity": {{
    "hepatotoxicity": "肝毒性风险",
    "cardiotoxicity": "心脏毒性风险",
    "mutagenicity": "致突变性"
  }}
}}
"""
        
        else:
            prompt = target
        
        return prompt
