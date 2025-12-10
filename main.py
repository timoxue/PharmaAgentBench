"""
PharmaAgentBench - 医药大模型评测系统
基于 AgentScope 多智能体框架
"""

try:
    import agentscope
    from agentscope.message import Msg
    AGENTSCOPE_AVAILABLE = True
except ImportError:
    print("警告: AgentScope 未正确安装，将使用基础模式")
    AGENTSCOPE_AVAILABLE = False
    
    # 使用 orchestrator_agent 中定义的基础类
    class Msg:
        def __init__(self, name: str, content: str, role: str = "assistant"):
            self.name = name
            self.content = content
            self.role = role

# 可选的 Web 框架支持
try:
    from flask import Flask, render_template_string, request
    FLASK_AVAILABLE = True
except ImportError:
    FLASK_AVAILABLE = False
from config import get_config
from agents.orchestrator_agent import OrchestratorAgent
import sys
import json
from datetime import datetime


def print_banner():
    """打印项目横幅"""
    banner = """
    ╔═══════════════════════════════════════════════════════════════╗
    ║                                                               ║
    ║           🧬 PharmaAgentBench 医药大模型评测系统 🧬            ║
    ║                                                               ║
    ║   基于 AgentScope 多智能体框架                                 ║
    ║   对比评测: DeepSeek-V3 | Qwen-Max | Gemini 3                 ║
    ║                                                               ║
    ╚═══════════════════════════════════════════════════════════════╝
    """
    print(banner)


def check_api_keys():
    """检查 API 密钥配置状态"""
    print("\n🔑 检查 API 密钥配置...")
    config = get_config()
    validation_results = config.validate_api_keys()
    
    all_valid = True
    for model, is_valid in validation_results.items():
        status = "✓" if is_valid else "✗"
        print(f"  {status} {model.upper()}: {'已配置' if is_valid else '未配置'}")
        if not is_valid:
            all_valid = False
    
    if not all_valid:
        print("\n⚠️  警告: 部分 API 密钥未配置，请在 config/api_keys.yaml 中填写")
        print("   或通过环境变量设置: DEEPSEEK_API_KEY, DASHSCOPE_API_KEY, GEMINI_API_KEY")
        response = input("\n是否继续运行? (y/N): ")
        if response.lower() != 'y':
            sys.exit(0)
    
    return validation_results


def initialize_agentscope():
    """初始化 AgentScope"""
    print("\n🚀 初始化 AgentScope 框架...")
    config = get_config()
    
    if not AGENTSCOPE_AVAILABLE:
        print("⚠️  AgentScope 未安装，跳过初始化")
        return
    
    try:
        # AgentScope 1.0.9+ 的简化初始化
        agentscope.init(
            project="PharmaAgentBench",
            name="pharma_eval_run",
        )
        
        print("✓ AgentScope 初始化完成")
    except Exception as e:
        print(f"⚠️  AgentScope 初始化失败: {e}")
        print("继续使用基础模式...")


def run_evaluation_task(task_config: dict):
    """运行评测任务
    
    Args:
        task_config: 任务配置
            - task_type: 任务类型 (molecule_generation, literature_extraction, admet_prediction)
            - target: 目标描述 (如 "EGFR 抑制剂")
            - models: 要评测的模型列表 (默认全部)
            - output_format: 输出格式 (默认 html)
    """
    print(f"\n📋 任务配置:")
    print(f"  类型: {task_config.get('task_type', 'molecule_generation')}")
    print(f"  目标: {task_config.get('target', 'N/A')}")
    print(f"  模型: {', '.join(task_config.get('models', ['deepseek', 'qwen', 'gemini']))}")
    
    # 创建主控智能体
    print("\n🤖 启动多智能体系统...")
    orchestrator = OrchestratorAgent(
        name="Orchestrator",
        model_config_name="gpt-4",  # 用于协调的模型
    )
    
    # 构造输入消息
    user_msg = Msg(
        name="User",
        content=json.dumps(task_config, ensure_ascii=False),
        role="user"
    )
    
    # 执行任务
    print("\n⚙️  执行评测任务...\n")
    result = orchestrator(user_msg)
    
    # 显示结果
    print("\n" + "="*70)
    print("📊 评测结果:")
    print("="*70)
    print(result.content)
    
    return result


def run_demo():
    """运行演示示例"""
    print("\n🎯 运行演示任务: 质子泵抑制剂分子设计")
    
    demo_task = {
        "task_type": "molecule_generation",
        "target": "模拟某已上市质子泵抑制剂的发现场景：靶点：胃壁细胞 H+/K+-ATPase，用于治疗胃酸相关疾病。要求：1. 分子量 300-450，LogP 1.5-3.5；2. QED >= 0.5，口服生物利用度良好；3. 避免强反应性官能团，并给出 SMILES、理化性质和设计理由。",
        "models": ["deepseek", "qwen", "gemini"],
        "constraints": {
            "max_mw": 450,  # 最大分子量
            "min_qed": 0.5,  # 最小 QED 评分
            "check_toxicity": True
        },
        "output_format": "html"
    }
    
    result = run_evaluation_task(demo_task)
    
    # 保存报告
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/output_{timestamp}.html"
    print(f"\n💾 报告已保存至: {report_path}")


def run_interactive():
    """交互式运行模式"""
    print("\n💬 进入交互式模式 (输入 'quit' 退出)")
    
    while True:
        print("\n" + "="*70)
        print("请选择任务类型:")
        print("  1. 分子生成 (molecule_generation)")
        print("  2. 文献信息提取 (literature_extraction)")
        print("  3. ADMET 预测 (admet_prediction)")
        print("  4. 退出")
        
        choice = input("\n请输入选项 (1-4): ").strip()
        
        if choice == '4' or choice.lower() == 'quit':
            print("\n👋 再见!")
            break
        
        task_types = {
            '1': 'molecule_generation',
            '2': 'literature_extraction',
            '3': 'admet_prediction'
        }
        
        if choice not in task_types:
            print("❌ 无效选项，请重新输入")
            continue
        
        task_type = task_types[choice]
        target = input("请输入任务目标描述: ").strip()
        
        if not target:
            print("❌ 目标描述不能为空")
            continue
        
        task_config = {
            "task_type": task_type,
            "target": target,
            "models": ["deepseek", "qwen", "gemini"],
            "output_format": "html"
        }
        
        try:
            run_evaluation_task(task_config)
        except Exception as e:
            print(f"\n❌ 任务执行失败: {str(e)}")
            import traceback
            traceback.print_exc()


def run_web_app(host: str = "0.0.0.0", port: int = 8000):
    """运行 Web 评测应用"""
    if not FLASK_AVAILABLE:
        print("\n❌ 未安装 Flask，无法启动 Web 应用。")
        print("   请先安装: pip install flask")
        return
    
    app = Flask(__name__)
    
    index_tpl = """<!DOCTYPE html>
<html lang=\"zh-CN\">
<head>
  <meta charset=\"UTF-8\">
  <title>PharmaAgentBench 在线评测</title>
  <style>
    body { font-family: -apple-system,BlinkMacSystemFont,'Segoe UI','Microsoft YaHei',sans-serif; background:#f5f5f5; padding:40px; }
    .container { max-width:900px; margin:0 auto; background:#fff; padding:30px; border-radius:10px; box-shadow:0 4px 20px rgba(0,0,0,0.08); }
    h1 { margin-bottom:10px; }
    label { font-weight:bold; }
    textarea { width:100%; height:120px; padding:10px; font-size:14px; }
    button { margin-top:15px; padding:10px 20px; font-size:15px; background:#667eea; color:#fff; border:none; border-radius:4px; cursor:pointer; }
    button:hover { background:#5563c1; }
    .meta { margin-top:20px; padding:10px; background:#f0f4ff; border-radius:6px; }
    pre { background:#111; color:#0f0; padding:15px; border-radius:6px; max-height:260px; overflow:auto; font-size:12px; }
    .summary { margin-top:20px; padding:15px; background:#f8f9fa; border-radius:6px; border-left:4px solid #667eea; white-space:pre-wrap; }
    a.btn { display:inline-block; margin-top:10px; padding:8px 16px; background:#28a745; color:#fff; text-decoration:none; border-radius:4px; }
  </style>
</head>
<body>
  <div class=\"container\">
    <h1>🧬 PharmaAgentBench 在线评测</h1>
    <p>输入你的任务描述，点击开始评测，对比 DeepSeek / Qwen / Gemini 在该任务下的表现。</p>
    <form method=\"post\">
      <label for=\"target\">目标描述：</label><br/>
      <textarea id=\"target\" name=\"target\">{{ default_target }}</textarea>
      <br/>
      <button type=\"submit\">开始评测</button>
    </form>
    {% if running %}
    <div class=\"meta\">
      <strong>评测状态：</strong> 已完成<br/>
      <strong>评测时间：</strong> {{ timestamp }}
    </div>
    <div class=\"summary\">
      {{ summary }}
    </div>
    <p>
      <a class=\"btn\" href=\"{{ report_path }}\" target=\"_blank\">打开完整 HTML 报告</a>
    </p>
    <h3>执行日志（用于查看进度）</h3>
    <pre>{{ log_text }}</pre>
    {% endif %}
  </div>
</body>
</html>"""
    
    @app.route('/', methods=['GET', 'POST'])
    def index():
        from datetime import datetime
        import io
        import contextlib
        
        default_target = "模拟某已上市质子泵抑制剂的发现场景：靠点：胃壁细胞 H+/K+-ATPase，用于治疗胃酸相关疾病。要求：1. 分子量 300-450，LogP 1.5-3.5；2. QED >= 0.5，口服生物利用度良好；3. 避免强反应性官能团，并给出 SMILES、理化性质和设计理由。"
        
        if request.method == 'POST':
            target = request.form.get('target', '').strip()
            if not target:
                target = default_target
            
            task_config = {
                "task_type": "molecule_generation",
                "target": target,
                "models": ["deepseek", "qwen", "gemini"],
                "constraints": {
                    "max_mw": 500,
                    "min_qed": 0.5,
                    "check_toxicity": True
                },
                "output_format": "html"
            }
            
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                result_msg = run_evaluation_task(task_config)
            log_text = buf.getvalue()
            
            content = result_msg.content if hasattr(result_msg, 'content') else str(result_msg)
            report_path = None
            summary = content
            for line in content.splitlines():
                if '完整报告已保存至' in line:
                    report_path = line.split('完整报告已保存至:')[-1].strip()
                    break
            
            return render_template_string(
                index_tpl,
                default_target=target,
                running=True,
                timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                summary=summary,
                report_path=report_path or '#',
                log_text=log_text
            )
        
        return render_template_string(
            index_tpl,
            default_target=default_target,
            running=False,
            timestamp="",
            summary="",
            report_path="#",
            log_text=""
        )
    
    print(f"\n🌐 Web 评测应用已启动: http://{host}:{port}")
    app.run(host=host, port=port, debug=False)


def main():
    """主函数"""
    print_banner()
    
    # 检查 API 密钥
    api_status = check_api_keys()
    
    # 初始化 AgentScope
    try:
        initialize_agentscope()
    except Exception as e:
        print(f"\n❌ AgentScope 初始化失败: {str(e)}")
        print("   提示: 请先安装 agentscope: pip install agentscope")
        sys.exit(1)
    
    # 解析命令行参数
    if len(sys.argv) > 1:
        if sys.argv[1] == '--demo':
            run_demo()
        elif sys.argv[1] == '--help':
            print("\n使用方法:")
            print("  python main.py           # 交互式模式")
            print("  python main.py --demo    # 运行演示示例")
            print("  python main.py --web     # 启动 Web 在线评测应用")
            print("  python main.py --help    # 显示帮助信息")
        elif sys.argv[1] == '--web':
            run_web_app()
        else:
            print(f"\n❌ 未知参数: {sys.argv[1]}")
            print("   使用 --help 查看帮助信息")
    else:
        # 默认交互式模式
        run_interactive()


if __name__ == "__main__":
    main()
