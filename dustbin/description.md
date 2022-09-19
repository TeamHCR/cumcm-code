## 代码说明

### 伪代码

```python
def generate_new():
    return 根据温度高低生成不同随机范围的位置

def run_random():
    # 降温到目标温度之前
    while T > T_finish:
        f_now = map(func(p), p in positions)
        next_positions = generate_new()
        for i in range(len(next_positions)):
            f = f_now[i]
            f_new = func(p)
            if metrospolis(f, f_new):
                positions[i] = next_positions[i]
        # 不断降温冷却
        T *= alpha
    # 返回点集中最好的
    return best(positions)

def run_climb_random():
    # 记录当前最值
    last = 0
    while T > T_finish:
        next_positions = generate_new()
        for p in next_positions:
            f = func(p)
            if f better than last:
                last = f; position = p;
        # 不断降温冷却
        T *= alpha
    return position
```

本文中为了求不可导的函数的最值，主要考虑了两种方法：模拟退火算法以及经过随机化修饰的爬山法。模拟退火算法在求解这类只能求得离散数值解的函数最值问题方面有很大的优势，其能够比较全局地考虑整个函数值域上的最值，且能够较快地解决最值问题。

但是，在本文的问题情景中，使用模拟退火算法进行一次求函数值的操作需要的时间接近 2 秒，且并不能在足够少的求值区间内得到函数的极值。尽管已经尽可能对模拟退火算法进行并行优化，如多线程、多进程、多服务器运行，其速度仍然达不到求解问题的需求。

相比之下，爬山法虽然只能通过周边函数值的启发，提高搜索的效率，但是其不是全面搜索，有遗漏函数极值的可能。

观察本问题情景中所求函数经过计算机插值之后生成的图像，可以发现其在大趋势上的变化很简单，这给了爬山法应用的条件；但是放大其图像细节则能够发现许多函数值的“山峰”，这些“山峰”很有可能阻挡爬山法向函数极值前进的步伐，但是模拟退火算法因其模拟的物理退火性质，能够随机地“越过”这些“山峰”，从而能够更加全面地发现函数的极值。

由模拟退火算法、爬山法和随机梯度下降算法等算法的启发，本文综合运用了以上多种算法，将由周围函数值启发前进方向、并行运行优化、变化的“温度”和随机跃进等思想结合起来，较好地完成了本文情景下的任务，不仅性能较为优异，其计算时间也大大缩短。