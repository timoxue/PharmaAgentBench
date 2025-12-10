"""
相似性检查器 - 用于比较不同模型生成的分子
"""

from typing import List, Dict, Any, Optional


class SimilarityChecker:
    """相似性检查器 - 比较分子结构的相似度"""
    
    def __init__(self):
        """初始化相似性检查器"""
        self.rdkit_available = False
        
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, DataStructs
            
            self.Chem = Chem
            self.AllChem = AllChem
            self.DataStructs = DataStructs
            self.rdkit_available = True
            
        except ImportError:
            print("警告: RDKit 未安装，相似度计算功能将受限")
    
    def calculate_pairwise_similarity(
        self,
        smiles_list: List[str],
        method: str = 'tanimoto'
    ) -> Dict[str, Any]:
        """计算一组分子的两两相似度
        
        Args:
            smiles_list: SMILES 列表
            method: 相似度计算方法
        
        Returns:
            相似度矩阵和统计信息
        """
        if not self.rdkit_available:
            return {
                'error': 'RDKit 未安装',
                'matrix': [],
                'statistics': {}
            }
        
        n = len(smiles_list)
        matrix = [[0.0] * n for _ in range(n)]
        
        # 解析所有 SMILES
        mols = []
        for smiles in smiles_list:
            mol = self.Chem.MolFromSmiles(smiles)
            mols.append(mol)
        
        # 生成指纹
        fps = []
        for mol in mols:
            if mol:
                fp = self.AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append(fp)
            else:
                fps.append(None)
        
        # 计算相似度矩阵
        similarities = []
        
        for i in range(n):
            for j in range(i + 1, n):
                if fps[i] and fps[j]:
                    if method == 'tanimoto':
                        sim = self.DataStructs.TanimotoSimilarity(fps[i], fps[j])
                    elif method == 'dice':
                        sim = self.DataStructs.DiceSimilarity(fps[i], fps[j])
                    else:
                        sim = self.DataStructs.TanimotoSimilarity(fps[i], fps[j])
                    
                    matrix[i][j] = sim
                    matrix[j][i] = sim
                    similarities.append(sim)
        
        # 对角线设为 1.0
        for i in range(n):
            matrix[i][i] = 1.0
        
        # 统计信息
        statistics = {}
        if similarities:
            statistics = {
                'mean': sum(similarities) / len(similarities),
                'min': min(similarities),
                'max': max(similarities),
                'count': len(similarities)
            }
        
        return {
            'matrix': matrix,
            'statistics': statistics,
            'method': method
        }
    
    def find_most_similar(
        self,
        query_smiles: str,
        candidate_smiles: List[str],
        top_k: int = 5
    ) -> List[Dict[str, Any]]:
        """找到与查询分子最相似的候选分子
        
        Args:
            query_smiles: 查询分子的 SMILES
            candidate_smiles: 候选分子的 SMILES 列表
            top_k: 返回前 k 个最相似的
        
        Returns:
            最相似分子列表
        """
        if not self.rdkit_available:
            return []
        
        query_mol = self.Chem.MolFromSmiles(query_smiles)
        if not query_mol:
            return []
        
        query_fp = self.AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
        
        results = []
        
        for smiles in candidate_smiles:
            mol = self.Chem.MolFromSmiles(smiles)
            if mol:
                fp = self.AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                similarity = self.DataStructs.TanimotoSimilarity(query_fp, fp)
                
                results.append({
                    'smiles': smiles,
                    'similarity': similarity
                })
        
        # 按相似度排序
        results.sort(key=lambda x: x['similarity'], reverse=True)
        
        return results[:top_k]
    
    def calculate_diversity(
        self,
        smiles_list: List[str]
    ) -> Dict[str, Any]:
        """计算一组分子的多样性
        
        Args:
            smiles_list: SMILES 列表
        
        Returns:
            多样性指标
        """
        # 计算两两相似度
        sim_result = self.calculate_pairwise_similarity(smiles_list)
        
        if 'error' in sim_result:
            return sim_result
        
        stats = sim_result.get('statistics', {})
        
        # 多样性 = 1 - 平均相似度
        mean_similarity = stats.get('mean', 0)
        diversity = 1.0 - mean_similarity
        
        return {
            'diversity_score': diversity,
            'mean_similarity': mean_similarity,
            'min_similarity': stats.get('min', 0),
            'max_similarity': stats.get('max', 0),
            'interpretation': self._interpret_diversity(diversity)
        }
    
    def _interpret_diversity(self, diversity_score: float) -> str:
        """解释多样性评分
        
        Args:
            diversity_score: 多样性评分
        
        Returns:
            解释文本
        """
        if diversity_score > 0.7:
            return "高度多样化 - 分子结构差异很大"
        elif diversity_score > 0.5:
            return "中等多样化 - 分子结构有一定差异"
        elif diversity_score > 0.3:
            return "低度多样化 - 分子结构较为相似"
        else:
            return "极低多样化 - 分子结构非常相似"
    
    def compare_model_outputs(
        self,
        model_results: Dict[str, str]
    ) -> Dict[str, Any]:
        """比较不同模型生成的分子
        
        Args:
            model_results: 模型名称到 SMILES 的映射
        
        Returns:
            比较结果
        """
        model_names = list(model_results.keys())
        smiles_list = [model_results[name] for name in model_names]
        
        # 计算相似度矩阵
        sim_result = self.calculate_pairwise_similarity(smiles_list)
        
        if 'error' in sim_result:
            return sim_result
        
        matrix = sim_result['matrix']
        n = len(model_names)
        
        # 构造模型对比表
        comparisons = []
        
        for i in range(n):
            for j in range(i + 1, n):
                comparisons.append({
                    'model1': model_names[i],
                    'model2': model_names[j],
                    'similarity': matrix[i][j],
                    'interpretation': self._interpret_similarity(matrix[i][j])
                })
        
        return {
            'comparisons': comparisons,
            'diversity': self.calculate_diversity(smiles_list),
            'matrix': matrix,
            'model_names': model_names
        }
    
    def _interpret_similarity(self, similarity: float) -> str:
        """解释相似度评分
        
        Args:
            similarity: 相似度评分
        
        Returns:
            解释文本
        """
        if similarity > 0.85:
            return "非常相似 - 可能是相同或高度相似的分子"
        elif similarity > 0.7:
            return "高度相似 - 核心骨架相同"
        elif similarity > 0.5:
            return "中等相似 - 部分结构特征相同"
        elif similarity > 0.3:
            return "低度相似 - 仅少量结构特征相同"
        else:
            return "不相似 - 结构差异很大"
