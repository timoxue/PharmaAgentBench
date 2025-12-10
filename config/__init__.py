"""配置管理模块"""

import os
import yaml
from pathlib import Path
from typing import Dict, Any


class ConfigManager:
    """API 密钥和配置管理器"""
    
    def __init__(self, config_path: str = None):
        """初始化配置管理器
        
        Args:
            config_path: 配置文件路径，默认为 config/api_keys.yaml
        """
        if config_path is None:
            config_path = Path(__file__).parent / "api_keys.yaml"
        
        self.config_path = Path(config_path)
        self.config = self._load_config()
    
    def _load_config(self) -> Dict[str, Any]:
        """加载配置文件"""
        if not self.config_path.exists():
            raise FileNotFoundError(f"配置文件不存在: {self.config_path}")
        
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        # 支持从环境变量覆盖
        config = self._override_from_env(config)
        
        return config
    
    def _override_from_env(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """从环境变量覆盖配置"""
        # DeepSeek
        if os.getenv('DEEPSEEK_API_KEY'):
            config['deepseek']['api_key'] = os.getenv('DEEPSEEK_API_KEY')
        
        # Qwen
        if os.getenv('DASHSCOPE_API_KEY'):
            config['qwen']['api_key'] = os.getenv('DASHSCOPE_API_KEY')
        
        # Gemini
        if os.getenv('GEMINI_API_KEY'):
            config['gemini']['api_key'] = os.getenv('GEMINI_API_KEY')
        if os.getenv('GCP_PROJECT_ID'):
            config['gemini']['project_id'] = os.getenv('GCP_PROJECT_ID')
        
        return config
    
    def get(self, key: str, default=None) -> Any:
        """获取配置值
        
        Args:
            key: 配置键，支持点分隔符，如 'deepseek.api_key'
            default: 默认值
        
        Returns:
            配置值
        """
        keys = key.split('.')
        value = self.config
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        
        return value
    
    def get_deepseek_config(self) -> Dict[str, Any]:
        """获取 DeepSeek 配置"""
        return self.config.get('deepseek', {})
    
    def get_qwen_config(self) -> Dict[str, Any]:
        """获取 Qwen 配置"""
        return self.config.get('qwen', {})
    
    def get_gemini_config(self) -> Dict[str, Any]:
        """获取 Gemini 配置"""
        return self.config.get('gemini', {})
    
    def validate_api_keys(self) -> Dict[str, bool]:
        """验证所有 API 密钥是否已配置
        
        Returns:
            各模型的 API 密钥验证状态
        """
        results = {}
        
        # 检查 DeepSeek
        deepseek_key = self.get('deepseek.api_key', '')
        results['deepseek'] = bool(deepseek_key and deepseek_key != 'YOUR_DEEPSEEK_API_KEY')
        
        # 检查 Qwen
        qwen_key = self.get('qwen.api_key', '')
        results['qwen'] = bool(qwen_key and qwen_key != 'YOUR_QWEN_API_KEY')
        
        # 检查 Gemini
        gemini_key = self.get('gemini.api_key', '')
        results['gemini'] = bool(gemini_key and gemini_key != 'YOUR_GEMINI_API_KEY')
        
        return results


# 全局配置实例
_config = None


def get_config() -> ConfigManager:
    """获取全局配置实例"""
    global _config
    if _config is None:
        _config = ConfigManager()
    return _config
