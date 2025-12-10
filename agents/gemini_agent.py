"""
Gemini Agent - Google Gemini 模型智能体
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
import os
import time
from config import get_config


class GeminiAgent(AgentBase):
    """Gemini 模型智能体"""
    
    def __init__(
        self,
        name: str = "Gemini",
        **kwargs
    ):
        """初始化 Gemini 智能体
        
        Args:
            name: 智能体名称
            **kwargs: 其他参数
        """
        super().__init__(name=name, **kwargs)
        
        # 获取配置
        config = get_config()
        self.gemini_config = config.get_gemini_config()
        
        # 初始化 Gemini 客户端
        self.client = self._init_client()
        self.model_name = self.gemini_config.get('model', 'gemini-3-pro-preview')
        
        # 加载 Prompt 模板
        self.prompts = self._load_prompts()
    
    def _init_client(self):
        """初始化 Gemini 客户端"""
        try:
            from google import genai
            import ssl
            
            # 配置 API 密钥
            api_key = self.gemini_config.get('api_key')
            
            # 设置环境变量（新 SDK 要求）
            os.environ['GOOGLE_API_KEY'] = api_key
            
            # 创建客户端，增加 SSL 容错配置
            try:
                client = genai.Client(api_key=api_key)
            except Exception as e:
                # 如果 SSL 错误，尝试使用更宽松的配置
                print(f"  ⚠️  Gemini 客户端初始化警告: {e}")
                client = genai.Client(api_key=api_key)
            
            return client
            
        except ImportError:
            # 如果未安装 google-genai，使用 REST API 替代方案
            print("警告: 未安装 google-genai，将使用 REST API")
            print("安装方法: pip install google-genai")
            return None
    
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
        
        # 调用 Gemini API
        print(f"\n  🔍 [Gemini] 开始调用 API...")
        print(f"  🔍 [Gemini] 任务类型: {task_type}")
        print(f"  🔍 [Gemini] 目标描述: {target[:80]}..." if len(target) > 80 else f"  🔍 [Gemini] 目标描述: {target}")
        print(f"  🔍 [Gemini] 模型名称: {self.model_name}")
        
        try:
            # 优先使用 REST API，避免 SDK 的 SSL 问题
            print(f"  🔍 [Gemini] 尝试使用 REST API...")
            try:
                content = self._call_with_rest_api(system_prompt, user_prompt)
                print(f"  ✅ [Gemini] REST API 调用成功")
            except Exception as rest_error:
                print(f"  ⚠️  [Gemini] REST API 失败: {rest_error}")
                # 如果 REST API 失败，尝试 SDK
                if self.client is not None:
                    print(f"  🔄 [Gemini] 切换到 SDK 调用...")
                    content = self._call_with_sdk(system_prompt, user_prompt)
                    print(f"  ✅ [Gemini] SDK 调用成功")
                else:
                    print(f"  ❌ [Gemini] SDK 未初始化，无法重试")
                    raise rest_error
            
            print(f"  🔍 [Gemini] 返回内容长度: {len(content)} 字符")
            print(f"  🔍 [Gemini] 内容预览: {content[:200]}..." if len(content) > 200 else f"  🔍 [Gemini] 完整内容: {content}")
            
            # 封装结果
            result = {
                'model': 'Gemini-3-Pro',
                'task_type': task_type,
                'target': target,
                'output': content,
                'metadata': {}
            }
            
            return Msg(
                name=self.name,
                content=json.dumps(result, ensure_ascii=False),
                role="assistant"
            )
            
        except Exception as e:
            print(f"  ❌ [Gemini] 调用失败: {type(e).__name__}: {str(e)}")
            error_result = {
                'model': 'Gemini-3-Pro',
                'task_type': task_type,
                'error': str(e)
            }
            
            return Msg(
                name=self.name,
                content=json.dumps(error_result, ensure_ascii=False),
                role="assistant"
            )
    
    def _call_with_sdk(self, system_prompt: str, user_prompt: str) -> str:
        """使用 SDK 调用 Gemini
        
        Args:
            system_prompt: 系统提示
            user_prompt: 用户提示
        
        Returns:
            模型响应
        """
        # 合并系统提示和用户提示
        combined_prompt = f"{system_prompt}\n\n{user_prompt}"
        
        # 重试配置
        max_retries = 3
        base_delay = 2  # 基础延迟（秒）
        
        for attempt in range(max_retries):
            try:
                # 使用新的 SDK API
                response = self.client.models.generate_content(
                    model=self.model_name,
                    contents=combined_prompt
                )
                
                return response.text
                
            except Exception as e:
                error_str = str(e)
                error_type = type(e).__name__
                
                # SSL 错误处理
                if 'SSL' in error_str or 'ssl' in error_str or 'UNEXPECTED_EOF' in error_str:
                    if attempt < max_retries - 1:
                        delay = base_delay * (2 ** attempt)
                        print(f"  ⚠️  Gemini SSL 连接错误，{delay}秒后重试... (第 {attempt + 1}/{max_retries} 次)")
                        time.sleep(delay)
                        # 重新初始化客户端
                        try:
                            self.client = self._init_client()
                        except:
                            pass
                        continue
                    else:
                        print(f"  ❌ Gemini SSL 错误，已重试 {max_retries} 次，切换到 REST API")
                        # 尝试使用 REST API 作为备用
                        try:
                            return self._call_with_rest_api(system_prompt, user_prompt)
                        except:
                            raise e
                
                # 检查是否是服务器过载或限额错误
                elif '503' in error_str or 'UNAVAILABLE' in error_str or 'overloaded' in error_str:
                    if attempt < max_retries - 1:
                        # 指数退避：2秒, 4秒, 8秒
                        delay = base_delay * (2 ** attempt)
                        print(f"  ⚠️  Gemini 服务器过载，{delay}秒后重试... (第 {attempt + 1}/{max_retries} 次)")
                        time.sleep(delay)
                        continue
                    else:
                        print(f"  ❌ Gemini 重试 {max_retries} 次后仍然失败")
                        raise
                
                # 配额限制错误
                elif '429' in error_str or 'RESOURCE_EXHAUSTED' in error_str:
                    # 从错误信息中提取建议的重试时间
                    import re
                    retry_match = re.search(r'retry in ([\d.]+)s', error_str)
                    if retry_match and attempt < max_retries - 1:
                        retry_seconds = float(retry_match.group(1))
                        print(f"  ⚠️  Gemini 配额超限，{retry_seconds:.1f}秒后重试...")
                        time.sleep(retry_seconds)
                        continue
                    else:
                        raise
                
                # 网络连接错误
                elif 'Connection' in error_str or 'Timeout' in error_str or 'EOF' in error_str:
                    if attempt < max_retries - 1:
                        delay = base_delay * (2 ** attempt)
                        print(f"  ⚠️  Gemini 网络连接错误，{delay}秒后重试...")
                        time.sleep(delay)
                        continue
                    else:
                        # 最后一次尝试 REST API
                        try:
                            print(f"  🔄 尝试使用 REST API 作为备用...")
                            return self._call_with_rest_api(system_prompt, user_prompt)
                        except:
                            raise e
                else:
                    # 其他错误直接抛出
                    raise
        
        # 如果所有重试都失败
        raise Exception(f"Gemini API 调用失败，已重试 {max_retries} 次")
    
    def _call_with_rest_api(self, system_prompt: str, user_prompt: str) -> str:
        """使用 REST API 调用 Gemini (备用方案)
        
        Args:
            system_prompt: 系统提示
            user_prompt: 用户提示
        
        Returns:
            模型响应
        """
        import requests
        
        print(f"    🌐 [Gemini-REST] 准备 REST API 请求...")
        
        api_key = self.gemini_config.get('api_key')
        url = f"https://generativelanguage.googleapis.com/v1beta/models/{self.model_name}:generateContent?key={api_key}"
        
        print(f"    🌐 [Gemini-REST] API URL: {url[:80]}...")
        
        payload = {
            "contents": [{
                "parts": [{
                    "text": f"{system_prompt}\n\n{user_prompt}"
                }]
            }],
            "generationConfig": {
                "temperature": 0.7,
                "maxOutputTokens": 2000
            }
        }
        
        # 重试配置
        max_retries = 3
        base_delay = 2
        
        for attempt in range(max_retries):
            try:
                print(f"    💬 [Gemini-REST] 发送请求 (尝试 {attempt + 1}/{max_retries})...")
                response = requests.post(url, json=payload, timeout=30)
                response.raise_for_status()
                
                print(f"    ✅ [Gemini-REST] 收到响应，状态码: {response.status_code}")
                
                result = response.json()
                print(f"    🔍 [Gemini-REST] 解析响应 JSON...")
                
                text = result['candidates'][0]['content']['parts'][0]['text']
                print(f"    ✅ [Gemini-REST] 成功获取内容，长度: {len(text)}")
                return text
                
            except requests.exceptions.HTTPError as e:
                print(f"    ⚠️  [Gemini-REST] HTTP 错误: {e.response.status_code} - {e.response.text[:200]}")
                if e.response.status_code in [503, 429]:
                    if attempt < max_retries - 1:
                        delay = base_delay * (2 ** attempt)
                        print(f"    ⏳ [Gemini-REST] {delay}秒后重试...")
                        time.sleep(delay)
                        continue
                raise
            except Exception as e:
                print(f"    ⚠️  [Gemini-REST] 请求失败: {type(e).__name__}: {str(e)}")
                if attempt < max_retries - 1:
                    delay = base_delay * (2 ** attempt)
                    print(f"    ⏳ [Gemini-REST] {delay}秒后重试...")
                    time.sleep(delay)
                    continue
                raise
        
        raise Exception(f"Gemini REST API 调用失败，已重试 {max_retries} 次")
    
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
        return """You are Gemini, an advanced AI assistant specialized in pharmaceutical chemistry and drug discovery.

Your expertise includes:
1. Rational drug design and molecular optimization
2. ADMET property prediction
3. Structure-activity relationship (SAR) analysis
4. Computational chemistry and cheminformatics
5. Biomedical literature mining

Guidelines:
- Represent molecular structures in SMILES format
- Provide scientifically rigorous analysis
- Consider synthetic accessibility
- Assess drug-likeness and ADMET properties
- Return results in JSON format with clear structure

Deliver professional, accurate, and actionable insights."""
    
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
            prompt = f"""Molecular Design Task:

Objective: {target}

Constraints:
- Maximum molecular weight: {constraints.get('max_mw', 500)} Da
- Minimum QED score: {constraints.get('min_qed', 0.5)}
- Toxicity assessment: {'Required' if constraints.get('check_toxicity') else 'Not required'}

Please return a JSON response with:
{{
  "smiles": "SMILES representation",
  "name": "Proposed molecule name",
  "properties": {{
    "molecular_weight": value,
    "logp": value,
    "hbd": value,
    "hba": value,
    "tpsa": value,
    "rotatable_bonds": value
  }},
  "qed_score": value,
  "design_rationale": "Explanation of design strategy",
  "synthesis_route": "Brief synthetic route suggestion",
  "predicted_activity": "Expected biological activity"
}}
"""
        
        elif task_type == 'literature_extraction':
            prompt = f"""Literature Mining Task:

Target: {target}

Please extract and return JSON with:
{{
  "key_findings": ["Finding 1", "Finding 2"],
  "molecules": [
    {{"name": "Name", "smiles": "SMILES", "activity": "Description", "ic50": "value"}}
  ],
  "clinical_data": {{
    "phase": "Clinical phase",
    "indication": "Indication",
    "results": "Summary"
  }},
  "mechanisms": ["Mechanism 1", "Mechanism 2"],
  "references": ["DOI/PMID 1", "DOI/PMID 2"]
}}
"""
        
        elif task_type == 'admet_prediction':
            prompt = f"""ADMET Prediction Task:

Molecule: {target}

Please predict ADMET properties and return JSON:
{{
  "absorption": {{
    "caco2_permeability": "High/Medium/Low",
    "human_intestinal_absorption": "percentage",
    "bioavailability_score": value
  }},
  "distribution": {{
    "plasma_protein_binding": "percentage",
    "vd": "L/kg",
    "bbb_penetration": "Yes/No"
  }},
  "metabolism": {{
    "cyp_substrates": ["CYP3A4", "CYP2D6"],
    "cyp_inhibitors": ["List"],
    "metabolic_stability": "High/Medium/Low"
  }},
  "excretion": {{
    "clearance": "mL/min/kg",
    "half_life": "hours",
    "renal_excretion": "percentage"
  }},
  "toxicity": {{
    "hepatotoxicity": "Risk level",
    "cardiotoxicity": "hERG inhibition risk",
    "mutagenicity": "Ames test prediction",
    "ld50": "mg/kg"
  }},
  "overall_assessment": "Summary and recommendations"
}}
"""
        
        else:
            prompt = target
        
        return prompt
