"""
Reporter Agent - 报告生成智能体
负责生成可视化评测报告
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
from datetime import datetime
from pathlib import Path


class ReporterAgent(AgentBase):
    """报告生成智能体 - 生成可视化 HTML 报告"""
    
    def __init__(
        self,
        name: str = "Reporter",
        **kwargs
    ):
        """初始化报告生成智能体
        
        Args:
            name: 智能体名称
            **kwargs: 其他参数
        """
        super().__init__(name=name, **kwargs)
        
        self.reports_dir = Path(__file__).parent.parent / "reports"
        self.reports_dir.mkdir(exist_ok=True)
    
    def reply(self, x: Optional[Union[Msg, List[Msg]]] = None) -> Msg:
        """生成评测报告
        
        Args:
            x: 输入消息，包含模型结果和审核结果
        
        Returns:
            报告内容
        """
        if isinstance(x, list):
            x = x[-1]
        
        # 解析输入
        try:
            data = json.loads(x.content) if isinstance(x.content, str) else x.content
        except json.JSONDecodeError:
            return Msg(
                name=self.name,
                content="报告生成失败: 无效的输入格式",
                role="assistant"
            )
        
        model_results = data.get('model_results', {})
        audit_results = data.get('audit_results', {})
        task_config = data.get('task_config', {})
        
        # 生成 HTML 报告
        html_content = self._generate_html_report(
            model_results,
            audit_results,
            task_config
        )
        
        # 保存报告
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        task_type = task_config.get('task_type', 'unknown')
        report_filename = f"report_{task_type}_{timestamp}.html"
        report_path = self.reports_dir / report_filename
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        # 生成文本摘要
        summary = self._generate_text_summary(
            model_results,
            audit_results,
            task_config
        )
        
        result = f"""{summary}

📄 完整报告已保存至: {report_path}
"""
        
        return Msg(
            name=self.name,
            content=result,
            role="assistant"
        )
    
    def _generate_text_summary(
        self,
        model_results: Dict[str, Any],
        audit_results: Dict[str, Any],
        task_config: Dict[str, Any]
    ) -> str:
        """生成文本摘要
        
        Args:
            model_results: 模型结果
            audit_results: 审核结果
            task_config: 任务配置
        
        Returns:
            文本摘要
        """
        task_type = task_config.get('task_type', 'unknown')
        target = task_config.get('target', 'N/A')
        
        summary_parts = [
            "=" * 70,
            "📊 评测报告摘要",
            "=" * 70,
            f"\n【任务类型】{task_type}",
            f"【目标描述】{target}",
            f"【评测时间】{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "\n【模型表现】"
        ]
        
        # 统计各模型的状态
        for model_name, result in model_results.items():
            if result.get('success'):
                audit = audit_results.get('audit_results', {}).get(model_name, {})
                status = "✓ 通过验证" if audit.get('valid') else "✗ 未通过验证"
                warnings = len(audit.get('warnings', []))
                warning_text = f" ({warnings} 个警告)" if warnings > 0 else ""
                
                summary_parts.append(
                    f"  • {model_name.upper():12s}: {status}{warning_text}"
                )
            else:
                summary_parts.append(
                    f"  • {model_name.upper():12s}: ✗ 执行失败"
                )
        
        # 添加审核摘要
        audit_summary = audit_results.get('summary', {})
        if audit_summary:
            summary_parts.append("\n【审核结果】")
            summary_parts.append(
                f"  通过率: {audit_summary.get('pass_rate', 0) * 100:.0f}% "
                f"({audit_summary.get('valid_models', 0)}/{audit_summary.get('total_models', 0)})"
            )
            
            recommendations = audit_summary.get('recommendations', [])
            if recommendations:
                summary_parts.append("\n【建议】")
                for rec in recommendations:
                    summary_parts.append(f"  {rec}")
        
        summary_parts.append("\n" + "=" * 70)
        
        return "\n".join(summary_parts)
    
    def _generate_html_report(
        self,
        model_results: Dict[str, Any],
        audit_results: Dict[str, Any],
        task_config: Dict[str, Any]
    ) -> str:
        """生成 HTML 报告
        
        Args:
            model_results: 模型结果
            audit_results: 审核结果
            task_config: 任务配置
        
        Returns:
            HTML 内容
        """
        task_type = task_config.get('task_type', 'unknown')
        target = task_config.get('target', 'N/A')
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        
        # 生成模型对比表格
        model_comparison_html = self._generate_model_comparison_table(
            model_results,
            audit_results
        )
        
        # 生成多维度对比分析（新增）
        dimensional_comparison_html = self._generate_dimensional_comparison(
            model_results,
            audit_results
        )
        
        # 生成详细结果
        detailed_results_html = self._generate_detailed_results(
            model_results,
            audit_results
        )
        
        # 完整 HTML
        html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PharmaAgentBench 评测报告</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'PingFang SC', 'Hiragino Sans GB', 
                         'Microsoft YaHei', sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            line-height: 1.6;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}
        
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        
        .header p {{
            opacity: 0.9;
            font-size: 1.1em;
        }}
        
        .section {{
            padding: 30px 40px;
        }}
        
        .section h2 {{
            color: #667eea;
            margin-bottom: 20px;
            font-size: 1.8em;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
        }}
        
        .info-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .info-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        
        .info-card h3 {{
            color: #667eea;
            margin-bottom: 10px;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        
        .info-card p {{
            font-size: 1.2em;
            color: #333;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        
        th {{
            background: #667eea;
            color: white;
            font-weight: 600;
        }}
        
        tr:hover {{
            background: #f5f5f5;
        }}
        
        .status-badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 0.85em;
            font-weight: 600;
        }}
        
        .status-success {{
            background: #d4edda;
            color: #155724;
        }}
        
        .status-warning {{
            background: #fff3cd;
            color: #856404;
        }}
        
        .status-error {{
            background: #f8d7da;
            color: #721c24;
        }}
        
        .model-result {{
            background: #f8f9fa;
            border-radius: 8px;
            padding: 20px;
            margin-bottom: 20px;
        }}
        
        .model-result h3 {{
            color: #667eea;
            margin-bottom: 15px;
        }}
        
        .result-content {{
            background: white;
            padding: 15px;
            border-radius: 6px;
            border: 1px solid #e0e0e0;
            font-family: 'Courier New', monospace;
            white-space: pre-wrap;
            word-wrap: break-word;
            max-height: 400px;
            overflow-y: auto;
        }}
        
        .footer {{
            background: #f8f9fa;
            padding: 20px;
            text-align: center;
            color: #666;
            border-top: 1px solid #e0e0e0;
        }}
        
        .warning-list {{
            list-style: none;
            padding-left: 0;
        }}
        
        .warning-list li {{
            background: #fff3cd;
            padding: 10px;
            margin: 5px 0;
            border-radius: 4px;
            border-left: 4px solid #ffc107;
        }}
        
        /* 新增：多维度对比样式 */
        .comparison-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        
        .dimension-card {{
            background: #f8f9fa;
            border-radius: 8px;
            padding: 20px;
            border-left: 4px solid #667eea;
        }}
        
        .dimension-card h3 {{
            color: #667eea;
            margin-bottom: 15px;
            font-size: 1.2em;
        }}
        
        .metric-row {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 10px;
            margin: 5px 0;
            background: white;
            border-radius: 4px;
        }}
        
        .metric-name {{
            font-weight: 600;
            color: #333;
        }}
        
        .metric-value {{
            font-family: 'Courier New', monospace;
            padding: 4px 10px;
            border-radius: 4px;
            font-weight: 600;
        }}
        
        .metric-best {{
            background: #d4edda;
            color: #155724;
        }}
        
        .metric-good {{
            background: #d1ecf1;
            color: #0c5460;
        }}
        
        .metric-warning {{
            background: #fff3cd;
            color: #856404;
        }}
        
        .metric-error {{
            background: #f8d7da;
            color: #721c24;
        }}
        
        .analysis-box {{
            background: #e7f3ff;
            border-left: 4px solid #2196F3;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }}
        
        .analysis-box h4 {{
            color: #2196F3;
            margin-bottom: 10px;
        }}
        
        .comparison-table {{
            width: 100%;
            margin: 20px 0;
            border-collapse: collapse;
        }}
        
        .comparison-table th {{
            background: #667eea;
            color: white;
            padding: 12px;
            text-align: center;
        }}
        
        .comparison-table td {{
            padding: 12px;
            text-align: center;
            border: 1px solid #ddd;
        }}
        
        .comparison-table tr:hover {{
            background: #f5f5f5;
        }}
        
        .score-bar {{
            height: 20px;
            background: #e0e0e0;
            border-radius: 10px;
            overflow: hidden;
            margin: 5px 0;
        }}
        
        .score-fill {{
            height: 100%;
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            transition: width 0.3s ease;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 PharmaAgentBench 评测报告</h1>
            <p>医药大模型对比评测系统</p>
        </div>
        
        <div class="section">
            <h2>📋 任务信息</h2>
            <div class="info-grid">
                <div class="info-card">
                    <h3>任务类型</h3>
                    <p>{task_type}</p>
                </div>
                <div class="info-card">
                    <h3>目标描述</h3>
                    <p>{target}</p>
                </div>
                <div class="info-card">
                    <h3>评测时间</h3>
                    <p>{timestamp}</p>
                </div>
                <div class="info-card">
                    <h3>参评模型</h3>
                    <p>{len(model_results)} 个</p>
                </div>
            </div>
        </div>
        
        <div class="section">
            <h2>📊 模型对比</h2>
            {model_comparison_html}
        </div>
        
        <div class="section">
            <h2>🎯 多维度对比分析</h2>
            {dimensional_comparison_html}
        </div>
        
        <div class="section">
            <h2>📝 详细结果</h2>
            {detailed_results_html}
        </div>
        
        <div class="footer">
            <p>Generated by PharmaAgentBench | {timestamp}</p>
            <p>Powered by AgentScope Multi-Agent Framework</p>
        </div>
    </div>
</body>
</html>
"""
        
        return html
    
    def _generate_model_comparison_table(
        self,
        model_results: Dict[str, Any],
        audit_results: Dict[str, Any]
    ) -> str:
        """生成模型对比表格"""
        rows = []
        
        for model_name, result in model_results.items():
            audit = audit_results.get('audit_results', {}).get(model_name, {})
            
            # 状态
            if not result.get('success'):
                status = '<span class="status-badge status-error">执行失败</span>'
            elif audit.get('valid'):
                status = '<span class="status-badge status-success">✓ 通过</span>'
            else:
                status = '<span class="status-badge status-error">✗ 未通过</span>'
            
            # 警告数
            warnings = len(audit.get('warnings', []))
            warning_badge = f'<span class="status-badge status-warning">{warnings}</span>' if warnings > 0 else '-'
            
            # 错误数
            errors = len(audit.get('errors', []))
            error_badge = f'<span class="status-badge status-error">{errors}</span>' if errors > 0 else '-'
            
            rows.append(f"""
                <tr>
                    <td><strong>{model_name.upper()}</strong></td>
                    <td>{status}</td>
                    <td>{warning_badge}</td>
                    <td>{error_badge}</td>
                </tr>
            """)
        
        table = f"""
        <table>
            <thead>
                <tr>
                    <th>模型</th>
                    <th>状态</th>
                    <th>警告</th>
                    <th>错误</th>
                </tr>
            </thead>
            <tbody>
                {''.join(rows)}
            </tbody>
        </table>
        """
        
        return table
    
    def _generate_detailed_results(
        self,
        model_results: Dict[str, Any],
        audit_results: Dict[str, Any]
    ) -> str:
        """生成详细结果"""
        sections = []
        
        for model_name, result in model_results.items():
            audit = audit_results.get('audit_results', {}).get(model_name, {})
            
            # 获取输出内容
            if result.get('success'):
                try:
                    response_data = result.get('response', '')
                    
                    # 尝试多种解析方式
                    if isinstance(response_data, str):
                        try:
                            response_data = json.loads(response_data)
                        except:
                            # 如果不是 JSON，直接使用原始字符串
                            output = response_data
                            response_data = None
                    
                    if response_data and isinstance(response_data, dict):
                        output = response_data.get('output', '')
                        
                        # 如果 output 为空，尝试其他字段
                        if not output:
                            # 尝试直接使用整个 response_data
                            output = json.dumps(response_data, ensure_ascii=False, indent=2)
                        elif isinstance(output, str):
                            # 尝试格式化 JSON
                            try:
                                output_obj = json.loads(output)
                                output = json.dumps(output_obj, ensure_ascii=False, indent=2)
                            except:
                                pass  # 保持原始字符串
                        else:
                            output = json.dumps(output, ensure_ascii=False, indent=2)
                    elif not output or output == '':
                        # 如果还是空的，显示原始结果
                        output = json.dumps(result, ensure_ascii=False, indent=2)
                    
                except Exception as e:
                    output = f"解析失败: {str(e)}\n\n原始数据:\n{json.dumps(result, ensure_ascii=False, indent=2)}"
            else:
                output = f"错误: {result.get('error', 'Unknown error')}"
            
            # HTML 转义避免显示问题
            import html
            output = html.escape(output)
            
            # 生成警告列表
            warnings_html = ""
            if audit.get('warnings'):
                warnings_items = ''.join(
                    f'<li>{html.escape(w)}</li>' for w in audit['warnings']
                )
                warnings_html = f'<ul class="warning-list">{warnings_items}</ul>'
            
            section = f"""
            <div class="model-result">
                <h3>{model_name.upper()}</h3>
                {warnings_html}
                <h4>模型输出:</h4>
                <div class="result-content">{output}</div>
            </div>
            """
            
            sections.append(section)
        
        return '\n'.join(sections)
    
    def _generate_dimensional_comparison(
        self,
        model_results: Dict[str, Any],
        audit_results: Dict[str, Any]
    ) -> str:
        """生成多维度对比分析
        
        Args:
            model_results: 模型结果
            audit_results: 审核结果
        
        Returns:
            HTML 内容
        """
        import json
        import re
        
        # 提取各模型的分子数据
        molecules_data = {}
        for model_name, result in model_results.items():
            if result.get('success'):
                try:
                    response_data = result.get('response', '')
                    
                    # 解析 response
                    if isinstance(response_data, str):
                        try:
                            response_data = json.loads(response_data)
                        except:
                            pass
                    
                    # 提取 output 字段
                    if isinstance(response_data, dict):
                        output = response_data.get('output', '')
                        
                        # 尝试从 output 中提取 JSON
                        if isinstance(output, str):
                            # 查找 JSON block (支持嵌套)
                            json_matches = re.finditer(r'```json\s*({.*?})\s*```', output, re.DOTALL)
                            for match in json_matches:
                                try:
                                    mol_data = json.loads(match.group(1))
                                    if 'smiles' in mol_data:
                                        molecules_data[model_name] = mol_data
                                        break
                                except:
                                    pass
                            
                            # 如果没找到，尝试直接查找 JSON 对象
                            if model_name not in molecules_data:
                                # 查找包含 smiles 的 JSON
                                json_pattern = r'\{[^{}]*"smiles"[^{}]*\}'
                                simple_match = re.search(json_pattern, output)
                                if simple_match:
                                    try:
                                        mol_data = json.loads(simple_match.group())
                                        molecules_data[model_name] = mol_data
                                    except:
                                        pass
                        elif isinstance(output, dict) and 'smiles' in output:
                            molecules_data[model_name] = output
                except Exception as e:
                    print(f"  ⚠️  提取 {model_name} 分子数据失败: {e}")
        
        # 获取审核结果
        audit_data = audit_results.get('audit_results', {})
        
        # 生成对比表格
        comparison_html = self._generate_properties_comparison_table(molecules_data, audit_data)
        
        # 生成维度分析
        dimensional_analysis = self._generate_dimensional_analysis(molecules_data, audit_data, model_results)
        
        # 生成优劣势分析
        strengths_weaknesses = self._generate_strengths_weaknesses(molecules_data, audit_data, model_results)
        
        html = f"""
        <div class="analysis-box">
            <h4>🔍 对比说明</h4>
            <p>以下分析从多个维度对比三个模型的输出质量，帮助你了解每个模型的优势和不足。</p>
        </div>
        
        <h3>📋 分子性质对比</h3>
        {comparison_html}
        
        <h3>🎯 维度评分</h3>
        {dimensional_analysis}
        
        <h3>⚖️ 优劣势分析</h3>
        {strengths_weaknesses}
        """
        
        return html
    
    def _generate_properties_comparison_table(
        self,
        molecules_data: Dict[str, Any],
        audit_data: Dict[str, Any]
    ) -> str:
        """生成分子性质对比表格"""
        if not molecules_data:
            return "<p>⚠️ 暂无可对比的分子数据</p>"
        
        # 提取性质
        properties_to_compare = [
            ('smiles', 'SMILES 结构', 'text'),
            ('molecular_weight', '分子量 (Da)', 'number', 500),
            ('logp', 'LogP', 'number', 5),
            ('hbd', '氢键供体', 'number', 5),
            ('hba', '氢键受体', 'number', 10),
            ('tpsa', 'TPSA (Å²)', 'number', 140),
            ('qed_score', 'QED 评分', 'score', 1.0)
        ]
        
        rows = []
        rows.append("""
        <table class="comparison-table">
            <thead>
                <tr>
                    <th>性质</th>
                    <th>DeepSeek</th>
                    <th>Qwen</th>
                    <th>Gemini</th>
                    <th>最优</th>
                </tr>
            </thead>
            <tbody>
        """)
        
        for prop_key, prop_name, prop_type, *threshold in properties_to_compare:
            row_data = {'property': prop_name, 'values': {}, 'best': None}
            
            # 提取每个模型的值
            for model_name in ['deepseek', 'qwen', 'gemini']:
                value = 'N/A'
                if model_name in molecules_data:
                    mol_data = molecules_data[model_name]
                    if prop_key == 'smiles':
                        value = mol_data.get('smiles', 'N/A')[:30] + '...'
                    elif prop_key == 'qed_score':
                        value = mol_data.get('qed_score', 'N/A')
                    else:
                        props = mol_data.get('properties', {})
                        # 属性别名映射，兼容不同模型的字段命名
                        alias_map = {
                            'molecular_weight': ['molecular_weight', 'MW', 'MolWt'],
                            'logp': ['logp', 'LogP'],
                            'hbd': ['hbd', 'HBD'],
                            'hba': ['hba', 'HBA'],
                            'tpsa': ['tpsa', 'TPSA']
                        }
                        candidate_keys = alias_map.get(prop_key, [prop_key, prop_key.replace('_', '')])
                        for k in candidate_keys:
                            if k in props:
                                value = props[k]
                                break
                
                row_data['values'][model_name] = value
            
            # 确定最优值
            if prop_type == 'number' or prop_type == 'score':
                numeric_values = {}
                for model_name, value in row_data['values'].items():
                    try:
                        if value != 'N/A':
                            numeric_values[model_name] = float(value)
                    except:
                        pass
                
                if numeric_values:
                    if prop_key in ['molecular_weight', 'logp', 'hbd', 'hba', 'tpsa']:
                        # 越小越好（但不能过小）
                        if threshold:
                            valid_values = {k: v for k, v in numeric_values.items() if v <= threshold[0]}
                            if valid_values:
                                row_data['best'] = min(valid_values, key=valid_values.get)
                    else:
                        # 越大越好
                        row_data['best'] = max(numeric_values, key=numeric_values.get)
            
            # 生成表格行
            cells = [f"<td><strong>{prop_name}</strong></td>"]
            for model_name in ['deepseek', 'qwen', 'gemini']:
                value = row_data['values'][model_name]
                css_class = 'metric-best' if model_name == row_data['best'] else ''
                
                # 处理显示格式
                if prop_type in ['number', 'score'] and value != 'N/A':
                    try:
                        value = f"{float(value):.2f}"
                    except:
                        pass
                
                cells.append(f'<td class="{css_class}">{value}</td>')
            
            # 添加最优标记
            best_label = row_data['best'].upper() if row_data['best'] else '-'
            cells.append(f'<td>{best_label}</td>')
            
            rows.append(f"<tr>{''.join(cells)}</tr>")
        
        rows.append("""
            </tbody>
        </table>
        """)
        
        return ''.join(rows)
    
    def _generate_dimensional_analysis(
        self,
        molecules_data: Dict[str, Any],
        audit_data: Dict[str, Any],
        model_results: Dict[str, Any]
    ) -> str:
        """生成维度评分分析"""
        dimensions = [
            ('调用结果', 'call_status'),
            ('化学合法性', 'chemistry'),
            ('类药性', 'druglikeness'),
            ('毒性风险', 'toxicity')
        ]
        
        html_parts = ['<div class="comparison-grid">']
        
        for dim_name, dim_key in dimensions:
            scores = {}
            for model_name in ['deepseek', 'qwen', 'gemini']:
                audit = audit_data.get(model_name, {})
                checks = audit.get('checks', {})
                score = None
                
                if dim_key == 'call_status':
                    result = model_results.get(model_name, {})
                    if result:
                        score = 100 if result.get('success') else 0
                elif dim_key == 'chemistry':
                    chem_check = checks.get('chemistry', {})
                    if chem_check:
                        score = 100 if chem_check.get('valid') else 0
                elif dim_key == 'druglikeness':
                    drug_check = checks.get('druglikeness', {})
                    if drug_check:
                        lipinski = drug_check.get('lipinski_violations', 5)
                        score = max(0, 100 - lipinski * 20)
                elif dim_key == 'toxicity':
                    tox_check = checks.get('toxicity', {})
                    if tox_check:
                        alerts = len(tox_check.get('alerts', []))
                        score = max(0, 100 - alerts * 25)
                
                scores[model_name] = score
            
            html_parts.append(f"""
            <div class="dimension-card">
                <h3>{dim_name}</h3>
            """)
            
            for model_name in ['deepseek', 'qwen', 'gemini']:
                score = scores.get(model_name, None)
                if score is None:
                    css_class = ''
                    score_text = '未评估'
                    width = 0
                else:
                    css_class = 'metric-best' if score == 100 else 'metric-good' if score >= 80 else 'metric-warning' if score >= 60 else 'metric-error'
                    score_text = f"{score}分"
                    width = score
                
                html_parts.append(f"""
                <div class="metric-row">
                    <span class="metric-name">{model_name.upper()}</span>
                    <span class="metric-value {css_class}">{score_text}</span>
                </div>
                <div class="score-bar">
                    <div class="score-fill" style="width: {width}%"></div>
                </div>
                """)
            
            html_parts.append('</div>')
        
        html_parts.append('</div>')
        return ''.join(html_parts)
    
    def _generate_strengths_weaknesses(
        self,
        molecules_data: Dict[str, Any],
        audit_data: Dict[str, Any],
        model_results: Dict[str, Any]
    ) -> str:
        """生成优劣势分析"""
        analysis = {}
        
        for model_name in ['deepseek', 'qwen', 'gemini']:
            strengths = []
            weaknesses = []
            
            # 检查成功状态
            result = model_results.get(model_name, {})
            if not result.get('success'):
                weaknesses.append('模型调用失败')
                analysis[model_name] = {'strengths': strengths, 'weaknesses': weaknesses}
                continue
            
            # 检查审核结果
            audit = audit_data.get(model_name, {})
            
            if audit.get('valid'):
                strengths.append('✅ 通过化学合法性验证')
            else:
                weaknesses.append('❌ 分子结构不合法')
            
            # 检查类药性
            checks = audit.get('checks', {})
            drug_check = checks.get('druglikeness', {})
            if drug_check:
                violations = drug_check.get('lipinski_violations', 5)
                if violations == 0:
                    strengths.append('✅ 完全符合 Lipinski 规则')
                elif violations <= 1:
                    strengths.append('✅ 类药性良好')
                else:
                    weaknesses.append(f'⚠️ 有 {violations} 项 Lipinski 违反')
            
            # 检查毒性
            tox_check = checks.get('toxicity', {})
            if tox_check:
                alerts = tox_check.get('alerts', [])
                if len(alerts) == 0:
                    strengths.append('✅ 无毒性警告')
                else:
                    weaknesses.append(f'⚠️ 发现 {len(alerts)} 个毒性警告')
            
            # 检查 QED 评分
            if model_name in molecules_data:
                qed = molecules_data[model_name].get('qed_score', 0)
                try:
                    qed_float = float(qed)
                    if qed_float >= 0.7:
                        strengths.append(f'✅ QED 评分高 ({qed_float:.2f})')
                    elif qed_float >= 0.5:
                        strengths.append(f'✅ QED 评分合格 ({qed_float:.2f})')
                    else:
                        weaknesses.append(f'⚠️ QED 评分偏低 ({qed_float:.2f})')
                except:
                    pass
            
            analysis[model_name] = {
                'strengths': strengths if strengths else ['⚠️ 暂无明显优势'],
                'weaknesses': weaknesses if weaknesses else ['✅ 暂无明显不足']
            }
        
        # 生成 HTML
        html_parts = ['<div class="comparison-grid">']
        
        for model_name in ['deepseek', 'qwen', 'gemini']:
            model_analysis = analysis[model_name]
            
            html_parts.append(f"""
            <div class="dimension-card">
                <h3>{model_name.upper()}</h3>
                <h4 style="color: #28a745; margin-top: 10px;">✨ 优势</h4>
                <ul>
            """)
            
            for strength in model_analysis['strengths']:
                html_parts.append(f'<li>{strength}</li>')
            
            html_parts.append("""
                </ul>
                <h4 style="color: #dc3545; margin-top: 15px;">⚠️ 不足</h4>
                <ul>
            """)
            
            for weakness in model_analysis['weaknesses']:
                html_parts.append(f'<li>{weakness}</li>')
            
            html_parts.append("""
                </ul>
            </div>
            """)
        
        html_parts.append('</div>')
        return ''.join(html_parts)
