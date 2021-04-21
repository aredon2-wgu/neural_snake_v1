"""Genetic Algorithmn Implementation
"""
#引入随机函数
import random

#定义基因遗传算法顶层接口
class GeneticAlgorithm(object):

    # 定义算法初始化
    def __init__(self, genetics):
        self.genetics = genetics #基因复制到个体
        pass
        
    # 定义算法run方法
    def run(self):
        population = self.genetics.initial()#原始种群
        while True:
            fits_pops = [(self.genetics.fitness(ch),  ch) for ch in population]  #对于原始种群中适应选择的个体 定义为适合种群
            if self.genetics.check_stop(fits_pops): break#当满足停止条件后退出 ，否则进行种群繁殖
            population = self.next(fits_pops)
            pass
        return population
        
    # 繁殖 精英主义思想
    def next(self, fits):
        parents_generator = self.genetics.parents(fits)# 从适应种群中选出双亲
        size = len(fits) #大小等于适应种群的数量
        nexts = []#初始化下一代列表为空
        while len(nexts) < size:#如果新的一代的数量小于适应种群的数量
            parents = next(parents_generator)#继续繁殖
            cross = random.random() < self.genetics.probability_crossover()#交叉限定
            children = self.genetics.crossover(parents) if cross else parents# 双亲基因进行交叉，通过则为子代否则为双亲
            for ch in children: #对子代中的每个个体
                mutate = random.random() < self.genetics.probability_mutation()#突变限定
                nexts.append(self.genetics.mutation(ch) if mutate else ch) # 适合个体添加到next列表中
                pass
            pass
        return nexts[0:size] # 返回下一代个体种群
    pass

class GeneticFunctions(object):

    #定义发生交叉率
    def probability_crossover(self):
        r"""returns rate of occur crossover(0.0-1.0)"""
        return 1.0
        
    #定义发生突变率
    def probability_mutation(self):
        r"""returns rate of occur mutation(0.0-1.0)"""
        return 0.0
        
    #初始种群的返回列表
    def initial(self):
        r"""returns list of initial population
        """
        return []
        
    #选择 返回适值
    def fitness(self, chromosome):
        r"""returns domain fitness value of chromosome
        """
        return len(chromosome)
        
    #判断是否达到停止条件  返回真 停止繁殖 适应种群（适值，染色体）
    def check_stop(self, fits_populations):
        r"""stop run if returns True
        - fits_populations: list of (fitness_value, chromosome)
        """
        return False
        
    #定义生成双亲
    def parents(self, fits_populations):
        r"""generator of selected parents
        """
        # 使用迭代器访问集合内已排序的元素
        gen = iter(sorted(fits_populations))
        while True:
            f1, ch1 = next(gen)
            f2, ch2 = next(gen)
            yield (ch1, ch2)#返回生成器中的双亲
            pass
        return
        
    #定义交叉
    def crossover(self, parents):
        r"""breed children
        """
        return parents
        
    #定义突变
    def mutation(self, chromosome):
        r"""mutate chromosome
        """
        return chromosome
    pass

if __name__ == "__main__":
    """
    #映射准备的匹配文本
    example: Mapped guess prepared Text
    """
    class GuessText(GeneticFunctions):
        #初始化 个体，目标文本，限制代数200，种群大小400，交叉概率0.9，后代变异概率0.2
        def __init__(self, target_text,
                     limit=200, size=400,
                     prob_crossover=0.9, prob_mutation=0.2):
            
            self.target = self.text2chromo(target_text)#将目标文本转换为ASCII码
            self.counter = 0
            self.limit = limit
            self.size = size
            self.prob_crossover = prob_crossover
            self.prob_mutation = prob_mutation
            pass
            
        # 基因遗传算法接口的实现
        # GeneticFunctions interface impls
       
        # 自身交叉概率
        def probability_crossover(self):
            return self.prob_crossover
            
        # 后代突变概率
        def probability_mutation(self):
            return self.prob_mutation
            
        # 初始指定数量的随机种群
        def initial(self):
            return [self.random_chromo() for j in range(self.size)]
            
        #定义适值函数 选择标准
        def fitness(self, chromo):
            # larger is better, matched == 0
            return -sum(abs(c - t) for c, t in zip(chromo, self.target))# 计算公式,返回对象中对应元素打包成一个元组 (基因,个体)
        
        # 判断是否达到停止条件
        def check_stop(self, fits_populations):
            self.counter += 1 #自身计数+1
            if self.counter % 10 == 0: #如果自身计数是整10数
                best_match = list(sorted(fits_populations))[-1][1] #倒序列表适应的种群
                fits = [f for f, ch in fits_populations]#从适应种群中依次取出适应个体
                best = max(fits) #最优适应
                worst = min(fits)#最劣适应
                ave = sum(fits) / len(fits) #求适应平均值
                # 打印输入 
                #[代数：%3d] score=(最优，平均，最劣)当前染色体值变异匹配文本
                print(
                    "[G %3d] score=(%4d, %4d, %4d): %r" %
                    (self.counter, best, ave, worst,
                     self.chromo2text(best_match)))
                pass
            return self.counter >= self.limit#返回大于适值函数的个体
            
        # 定义双亲
        def parents(self, fits_populations):
            while True:
                father = self.tournament(fits_populations)# 在适应种群中随机生成双亲
                mother = self.tournament(fits_populations)
                yield (father, mother)# 返回生成器中的双亲
                pass
            pass
            
        #交叉
        def crossover(self, parents):
            father, mother = parents
            index1 = random.randint(1, len(self.target) - 2)
            index2 = random.randint(1, len(self.target) - 2)
            if index1 > index2: index1, index2 = index2, index1 #保证index1<index2
            child1 = father[:index1] + mother[index1:index2] + father[index2:] #交叉
            child2 = mother[:index1] + father[index1:index2] + mother[index2:] 
            return (child1, child2) #返回交叉后的子代
            
        #突变
        def mutation(self, chromosome):
            index = random.randint(0, len(self.target) - 1) #从当前匹配文本长度中随机取一个整数
            vary = random.randint(-5, 5)#随机生成-5～5之间的数作为突变基因
            mutated = list(chromosome)#从容器中取出基因
            mutated[index] += vary# 突变
            return mutated

        # internals
        #竞争
        def tournament(self, fits_populations):      
            alicef, alice = self.select_random(fits_populations) #从适应种群中随机取出4个个体
            bobf, bob = self.select_random(fits_populations)
            return alice if alicef > bobf else bob  #如果alicef>bobf 返回alice 否则返回bob
            
        #随机选择
        def select_random(self, fits_populations):
            return fits_populations[random.randint(0, len(fits_populations)-1)]
            
        #定义文本到染色体
        #参数是一个ascii字符，返回值是对应的十进制整数
        def text2chromo(self, text):
            return [ord(ch) for ch in text]
            
        #定义染色体到文本
        #参数是0 - 256 的一个整数，返回值是当前整数对应的ascii字符。参数可以是10进制也可以是16进制的形式
        def chromo2text(self, chromo):
            return "".join(chr(max(1, min(ch, 255))) for ch in chromo)
            
        #定义随机
        def random_chromo(self):
            return [random.randint(1, 255) for i in range(len(self.target))]
        pass
        
    #利用基因遗传算法中进行逐渐匹配文本
    GeneticAlgorithm(GuessText("hello,world")).run()
    pass
