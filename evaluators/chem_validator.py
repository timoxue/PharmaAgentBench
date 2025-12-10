"""
化学验证器 - 使用 RDKit 验证分子结构
"""

from typing import Dict, Any, List, Optional


class ChemValidator:
    """化学验证器 - 验证分子合法性、类药性和毒性"""
    
    def __init__(self):
        """初始化化学验证器"""
        self.rdkit_available = False
        self.mol_descriptors = None
        self.qed = None
        
        # 尝试导入 RDKit
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, QED, Lipinski
            from rdkit.Chem import AllChem
            
            self.Chem = Chem
            self.Descriptors = Descriptors
            self.QED = QED
            self.Lipinski = Lipinski
            self.AllChem = AllChem
            self.rdkit_available = True
            
        except ImportError:
            print("警告: RDKit 未安装，化学验证功能将受限")
            print("安装方法: pip install rdkit")
    
    def validate_smiles(self, smiles: str) -> Dict[str, Any]:
        """验证 SMILES 字符串的合法性
        
        Args:
            smiles: SMILES 字符串
        
        Returns:
            验证结果
        """
        result = {
            'valid': False,
            'smiles': smiles,
            'canonical_smiles': None,
            'errors': []
        }
        
        if not self.rdkit_available:
            result['errors'].append("RDKit 未安装")
            return result
        
        if not smiles or not isinstance(smiles, str):
            result['errors'].append("SMILES 为空或格式错误")
            return result
        
        # 尝试解析 SMILES
        try:
            mol = self.Chem.MolFromSmiles(smiles)
            
            if mol is None:
                result['errors'].append("SMILES 解析失败，可能包含无效字符")
                return result
            
            # 生成规范 SMILES
            canonical_smiles = self.Chem.MolToSmiles(mol)
            
            result['valid'] = True
            result['canonical_smiles'] = canonical_smiles
            result['num_atoms'] = mol.GetNumAtoms()
            result['num_bonds'] = mol.GetNumBonds()
            
        except Exception as e:
            result['errors'].append(f"验证过程出错: {str(e)}")
        
        return result
    
    def check_druglikeness(self, smiles: str) -> Dict[str, Any]:
        """检查分子的类药性
        
        Args:
            smiles: SMILES 字符串
        
        Returns:
            类药性评估结果
        """
        result = {
            'valid': False,
            'properties': {},
            'lipinski_violations': 0,
            'veber_violations': 0,
            'qed_score': 0.0,
            'warnings': []
        }
        
        if not self.rdkit_available:
            result['warnings'].append("RDKit 未安装")
            return result
        
        try:
            mol = self.Chem.MolFromSmiles(smiles)
            if mol is None:
                result['warnings'].append("无效的 SMILES")
                return result
            
            # 计算分子性质
            properties = {
                'molecular_weight': self.Descriptors.MolWt(mol),
                'logp': self.Descriptors.MolLogP(mol),
                'hbd': self.Descriptors.NumHDonors(mol),
                'hba': self.Descriptors.NumHAcceptors(mol),
                'tpsa': self.Descriptors.TPSA(mol),
                'rotatable_bonds': self.Descriptors.NumRotatableBonds(mol),
                'aromatic_rings': self.Lipinski.NumAromaticRings(mol),
                'num_atoms': mol.GetNumAtoms(),
            }
            
            result['properties'] = properties
            
            # Lipinski's Rule of Five
            lipinski_violations = 0
            if properties['molecular_weight'] > 500:
                lipinski_violations += 1
                result['warnings'].append("分子量 > 500 Da (Lipinski 违规)")
            
            if properties['logp'] > 5:
                lipinski_violations += 1
                result['warnings'].append("LogP > 5 (Lipinski 违规)")
            
            if properties['hbd'] > 5:
                lipinski_violations += 1
                result['warnings'].append("氢键供体 > 5 (Lipinski 违规)")
            
            if properties['hba'] > 10:
                lipinski_violations += 1
                result['warnings'].append("氢键受体 > 10 (Lipinski 违规)")
            
            result['lipinski_violations'] = lipinski_violations
            
            # Veber's Rules
            veber_violations = 0
            if properties['rotatable_bonds'] > 10:
                veber_violations += 1
                result['warnings'].append("旋转键 > 10 (Veber 违规)")
            
            if properties['tpsa'] > 140:
                veber_violations += 1
                result['warnings'].append("TPSA > 140 Ų (Veber 违规)")
            
            result['veber_violations'] = veber_violations
            
            # 计算 QED (Quantitative Estimate of Drug-likeness)
            try:
                qed_score = self.QED.qed(mol)
                result['qed_score'] = qed_score
            except:
                result['warnings'].append("QED 计算失败")
            
            result['valid'] = True
            
        except Exception as e:
            result['warnings'].append(f"类药性检查失败: {str(e)}")
        
        return result
    
    def check_toxicity(self, smiles: str) -> Dict[str, Any]:
        """检查分子的潜在毒性
        
        Args:
            smiles: SMILES 字符串
        
        Returns:
            毒性评估结果
        """
        result = {
            'valid': False,
            'alerts': [],
            'warnings': [],
            'pains_alerts': []
        }
        
        if not self.rdkit_available:
            result['warnings'].append("RDKit 未安装")
            return result
        
        try:
            mol = self.Chem.MolFromSmiles(smiles)
            if mol is None:
                result['warnings'].append("无效的 SMILES")
                return result
            
            # 检查 PAINS (Pan-Assay Interference Compounds)
            pains_patterns = self._get_pains_patterns()
            
            for pattern_name, pattern_smarts in pains_patterns.items():
                pattern = self.Chem.MolFromSmarts(pattern_smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    result['pains_alerts'].append(pattern_name)
            
            # 检查常见毒性基团
            toxic_patterns = self._get_toxic_patterns()
            
            for pattern_name, pattern_smarts in toxic_patterns.items():
                pattern = self.Chem.MolFromSmarts(pattern_smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    result['alerts'].append(pattern_name)
            
            result['valid'] = True
            
        except Exception as e:
            result['warnings'].append(f"毒性检查失败: {str(e)}")
        
        return result
    
    def _get_pains_patterns(self) -> Dict[str, str]:
        """获取 PAINS 模式 (简化版)
        
        Returns:
            PAINS 模式字典
        """
        # 这里仅包含部分常见 PAINS 模式
        return {
            'Quinone': 'C1=CC(=O)C=CC1=O',
            'Michael_acceptor': 'C=CC(=O)',
            'Catechol': 'c1c(O)c(O)ccc1',
            'Rhodanine': 'C1C(=O)NC(=S)S1',
            'Alkyl_halide': '[CX4][F,Cl,Br,I]',
        }
    
    def _get_toxic_patterns(self) -> Dict[str, str]:
        """获取常见毒性基团模式
        
        Returns:
            毒性模式字典
        """
        return {
            'Nitro_aromatic': '[cR]N(=O)=O',
            'Aromatic_amine': 'cN',
            'Aldehyde': '[CH]=O',
            'Hydrazine': 'NN',
            'Azide': 'N=[N+]=[N-]',
            'Peroxide': 'OO',
            'Epoxide': 'C1OC1',
        }
    
    def calculate_similarity(
        self,
        smiles1: str,
        smiles2: str,
        method: str = 'tanimoto'
    ) -> float:
        """计算两个分子的相似度
        
        Args:
            smiles1: 第一个 SMILES
            smiles2: 第二个 SMILES
            method: 相似度计算方法 (tanimoto, dice)
        
        Returns:
            相似度分数 (0-1)
        """
        if not self.rdkit_available:
            return 0.0
        
        try:
            mol1 = self.Chem.MolFromSmiles(smiles1)
            mol2 = self.Chem.MolFromSmiles(smiles2)
            
            if mol1 is None or mol2 is None:
                return 0.0
            
            # 生成分子指纹
            fp1 = self.AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
            fp2 = self.AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
            
            # 计算相似度
            from rdkit import DataStructs
            
            if method == 'tanimoto':
                similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            elif method == 'dice':
                similarity = DataStructs.DiceSimilarity(fp1, fp2)
            else:
                similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            
            return similarity
            
        except Exception as e:
            print(f"相似度计算失败: {e}")
            return 0.0
