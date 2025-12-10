"""
Orchestrator Agent - 主控智能体
负责任务分发、结果汇总和流程编排
"""

try:
    from agentscope.agents import AgentBase
    from agentscope.message import Msg
except ImportError:
    # 尝试旧版本的导入方式
    try:
        from agentscope import AgentBase
        from agentscope import Msg
    except ImportError:
        # 如果都失败，提供一个基础实现
        class AgentBase:
            def __init__(self, name: str = "Agent", **kwargs):
                self.name = name
            
            def reply(self, x):
                raise NotImplementedError("Subclass must implement reply method")
            
            def __call__(self, x):
                return self.reply(x)
            
            def speak(self, msg):
                print(f"[{self.name}] {msg.content if hasattr(msg, 'content') else msg}")
        
        class Msg:
            def __init__(self, name: str, content: str, role: str = "assistant"):
                self.name = name
                self.content = content
                self.role = role
from typing import Optional, Union, List, Dict, Any
import json
from config import get_config


class OrchestratorAgent(AgentBase):
    """主控智能体 - 负责任务协调和流程编排"""
    
    def __init__(
        self,
        name: str = "Orchestrator",
        model_config_name: str = None,
        **kwargs
    ):
        """初始化主控智能体
        
        Args:
            name: 智能体名称
            model_config_name: 模型配置名称（可选）
            **kwargs: 其他参数
        """
        super().__init__(
            name=name,
            model_config_name=model_config_name,
            **kwargs
        )
        
        self.config = get_config()
        self.model_agents = {}
        self.auditor = None
        self.reporter = None
        
        # 延迟初始化子智能体（避免循环导入）
        self._agents_initialized = False
    
    def _initialize_agents(self):
        """初始化所有子智能体"""
        if self._agents_initialized:
            return
        
        from .deepseek_agent import DeepSeekAgent
        from .qwen_agent import QwenAgent
        from .gemini_agent import GeminiAgent
        from .auditor_agent import AuditorAgent
        from .reporter_agent import ReporterAgent
        
        # 检查哪些模型的 API 已配置
        api_status = self.config.validate_api_keys()
        
        # 初始化模型智能体
        if api_status.get('deepseek'):
            self.model_agents['deepseek'] = DeepSeekAgent(name="DeepSeek")
        
        if api_status.get('qwen'):
            self.model_agents['qwen'] = QwenAgent(name="Qwen")
        
        if api_status.get('gemini'):
            self.model_agents['gemini'] = GeminiAgent(name="Gemini")
        
        # 初始化审核和报告智能体
        self.auditor = AuditorAgent(name="Auditor")
        self.reporter = ReporterAgent(name="Reporter")
        
        self._agents_initialized = True
        
        self.speak(Msg(
            name=self.name,
            content=f"✓ 已初始化 {len(self.model_agents)} 个模型智能体",
            role="assistant"
        ))
    
    def reply(self, x: Optional[Union[Msg, List[Msg]]] = None) -> Msg:
        """处理消息并协调各智能体
        
        Args:
            x: 输入消息
        
        Returns:
            处理结果消息
        """
        # 初始化子智能体
        self._initialize_agents()
        
        # 解析任务配置
        if isinstance(x, list):
            x = x[-1]
        
        try:
            task_config = json.loads(x.content) if isinstance(x.content, str) else x.content
        except json.JSONDecodeError:
            # 如果不是 JSON，作为自然语言处理
            task_config = {
                "task_type": "molecule_generation",
                "target": x.content,
                "models": list(self.model_agents.keys())
            }
        
        task_type = task_config.get('task_type', 'molecule_generation')
        target = task_config.get('target', '')
        selected_models = task_config.get('models', list(self.model_agents.keys()))
        
        self.speak(Msg(
            name=self.name,
            content=f"📋 开始执行任务: {task_type}\n   目标: {target}",
            role="assistant"
        ))
        
        # 步骤 1: 分发任务给各模型智能体
        model_results = self._distribute_to_models(task_config, selected_models)
        
        # 步骤 2: 审核结果
        audit_results = self._audit_results(model_results, task_config)
        
        # 步骤 3: 生成报告
        report = self._generate_report(model_results, audit_results, task_config)
        
        return Msg(
            name=self.name,
            content=report,
            role="assistant"
        )
    
    def _distribute_to_models(
        self, 
        task_config: Dict[str, Any],
        selected_models: List[str]
    ) -> Dict[str, Any]:
        """分发任务给各模型智能体
        
        Args:
            task_config: 任务配置
            selected_models: 选定的模型列表
        
        Returns:
            各模型的返回结果
        """
        self.speak(Msg(
            name=self.name,
            content=f"🔄 分发任务给 {len(selected_models)} 个模型...",
            role="assistant"
        ))
        
        results = {}
        
        for model_name in selected_models:
            if model_name not in self.model_agents:
                self.speak(Msg(
                    name=self.name,
                    content=f"⚠️  跳过未配置的模型: {model_name}",
                    role="assistant"
                ))
                continue
            
            agent = self.model_agents[model_name]
            
            # 构造任务消息
            task_msg = Msg(
                name=self.name,
                content=json.dumps(task_config, ensure_ascii=False),
                role="user"
            )
            
            # 调用模型智能体
            try:
                self.speak(Msg(
                    name=self.name,
                    content=f"  → 调用 {model_name.upper()}...",
                    role="assistant"
                ))
                
                response = agent(task_msg)
                results[model_name] = {
                    'success': True,
                    'response': response.content,
                    'agent': model_name
                }
                
                self.speak(Msg(
                    name=self.name,
                    content=f"  ✓ {model_name.upper()} 完成",
                    role="assistant"
                ))
                
            except Exception as e:
                results[model_name] = {
                    'success': False,
                    'error': str(e),
                    'agent': model_name
                }
                
                self.speak(Msg(
                    name=self.name,
                    content=f"  ✗ {model_name.upper()} 失败: {str(e)}",
                    role="assistant"
                ))
        
        return results
    
    def _audit_results(
        self, 
        model_results: Dict[str, Any],
        task_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """审核各模型的结果
        
        Args:
            model_results: 各模型的返回结果
            task_config: 任务配置
        
        Returns:
            审核结果
        """
        self.speak(Msg(
            name=self.name,
            content="🔍 开始审核结果...",
            role="assistant"
        ))
        
        # 构造审核消息
        audit_msg = Msg(
            name=self.name,
            content=json.dumps({
                'model_results': model_results,
                'task_config': task_config
            }, ensure_ascii=False),
            role="user"
        )
        
        # 调用审核智能体
        try:
            audit_response = self.auditor(audit_msg)
            audit_results = json.loads(audit_response.content) if isinstance(audit_response.content, str) else audit_response.content
            
            self.speak(Msg(
                name=self.name,
                content="✓ 审核完成",
                role="assistant"
            ))
            
            return audit_results
            
        except Exception as e:
            self.speak(Msg(
                name=self.name,
                content=f"⚠️  审核失败: {str(e)}",
                role="assistant"
            ))
            
            return {'success': False, 'error': str(e)}
    
    def _generate_report(
        self,
        model_results: Dict[str, Any],
        audit_results: Dict[str, Any],
        task_config: Dict[str, Any]
    ) -> str:
        """生成评测报告
        
        Args:
            model_results: 各模型的返回结果
            audit_results: 审核结果
            task_config: 任务配置
        
        Returns:
            报告内容
        """
        self.speak(Msg(
            name=self.name,
            content="📊 生成评测报告...",
            role="assistant"
        ))
        
        # 构造报告消息
        report_msg = Msg(
            name=self.name,
            content=json.dumps({
                'model_results': model_results,
                'audit_results': audit_results,
                'task_config': task_config
            }, ensure_ascii=False),
            role="user"
        )
        
        # 调用报告智能体
        try:
            report_response = self.reporter(report_msg)
            
            self.speak(Msg(
                name=self.name,
                content="✓ 报告生成完成",
                role="assistant"
            ))
            
            return report_response.content
            
        except Exception as e:
            self.speak(Msg(
                name=self.name,
                content=f"⚠️  报告生成失败: {str(e)}",
                role="assistant"
            ))
            
            return f"报告生成失败: {str(e)}"
