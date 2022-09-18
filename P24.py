import threading
import queue
import time
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wl, wlexpr
import matplotlib.pyplot as plt
import random
import math
from functools import reduce
import nest_asyncio
nest_asyncio.apply()
# from wolframclient.evaluation import parallel_evaluate

P21 = """
Clear["Global`*"];
M = 4866; m = 2433; k = 80000; c = 10000; omega = 2.2143; k1 = 1025 * 9.8 * Pi;
f = 4890; c0 = 167.8395; m0 = 1165.992;
start = 510; stop = 520; dim = 0.01;
P21[cm2_] := Sum[cm2 * ((x1'[i] - x2'[i])^2) * dim, {i, start, stop, dim}] /. 
  NDSolve[{m*x2''[t]==k*(x1[t]-x2[t])+cm2*(x1'[t]-x2'[t]), 
    f*Cos[omega * t]==k*(x1[t]-x2[t])+cm2*(x1'[t]-x2'[t])+k1*x1[t]+c0*x1'[t]+m0*x1''[t]+M*x1''[t], 
    x1[0]==x2[0]==x1'[0]==x2'[0]==0},
  {x1, x2}, {t, start, stop}]
"""

P22 = """
Clear["Global`*"];
M = 4866; m = 2433; k = 80000; c = 10000; omega = 2.2143; k1 = 1025 * 9.8 * Pi;
f = 4890; c0 = 167.8395; m0 = 1165.992;
start = 110; stop = start + 20 * 1 / omega; dim = 0.1;
cmv2[v_, ck_, ce_] := ck * Abs[v]^ce;
P22I[ck_, ce_] := NDSolve[{m*x2''[t]==k*(x1[t]-x2[t])+cmv2[x1'[t]-x2'[t], ck, ce]*(x1'[t]-x2'[t]), 
    f*Cos[omega * t]==k*(x1[t]-x2[t])+cmv2[x1'[t]-x2'[t], ck, ce]*(x1'[t]-x2'[t])+k1*x1[t]+c0*x1'[t]+m0*x1''[t]+M*x1''[t], 
    x1[0]==x2[0]==x1'[0]==x2'[0]==0},
  {x1, x2}, {t, start, stop}]
P22[ck_, ce_] := Sum[cmv2[x1'[i]-x2'[i], ck, ce+2] * dim, {i, start, stop, dim}] /. P22I[ck, ce]
"""

P4 = '''
Clear["Global`*"];
(* AG长度 *)
d = 1.40792;
(* 弹簧原长度 *)
l0 = 0.5; 
(* 振子质量 *)
m = 2433; 
(* 弹簧劲度系数 *)
k = 80000;
(* 小g *)
g = 9.8;
(* 弹簧旋转刚度 *)
krot = 250000;
(* 浮子质量 *)
M = 4866;
(* 垂荡兴波阻尼系数 *)
c0 = 528.5018;
(* 纵摇附加转动惯量 *)
I0 = 7142.493;
(* 纵摇激励力矩振幅 *)
L = 2140;
(* 入射波浪频率 *)
w = 1.9806;
(* 垂荡激励力振幅 *)
f = 1760;
(* 浮子转动惯量 *)
(* Ia = 300*4866/(7 \[Pi]+(Sqrt[41] \[Pi])/5)*(Integrate[1/0.8*2*Pi*(r+(17.6+8/75*(Sqrt[41]))/(7+(Sqrt[41])/5))*r^2,{r,-(17.6+8/75*(Sqrt[41]))/(7+(Sqrt[41])/5),-(17.6+8/75*(Sqrt[41]))/(7+(Sqrt[41])/5)+0.8}]+Integrate[2*Pi*r^2,{r,-(17.6+8/75*(Sqrt[41]))/(7+(Sqrt[41])/5)+0.8,-(17.6+8/75*(Sqrt[41]))/(7+(Sqrt[41])/5)+3.8}]+Integrate[2*Pi*(3.8-(17.6+8/75*(Sqrt[41]))/(7+(Sqrt[41])/5))^2+a^2,{a,0,0.5}]); *)
Ia = 8399.2;

(* 垂荡附加质量 *)
m0 = 1091.099;
(* 纵摇兴波阻尼力矩系数 *)
cr0 = 1655.909;
(* 纵摇静水恢复力矩系数 *)
Mrec = 8890.7;
(* 海水密度 *)
rho = 1025;

(* 垂荡静水恢复力(浮力) *)
Fj[t_] := (7.120975609756098 - Pi*zg[t]) * rho * g;
(* 激励力 *)
Fwave[t_] := f * Cos[w*t];
(* 振子转动惯量 *)
(* Ib[t_] := 258.149 + 774.448 l[t] + 1548.9 l[t]^2; *)
(* Ib[t_] := MomentOfInertia[Cylinder[{{0, 0, 0}, {0.5, 0, 0}}, 0.5], {l[t] + 0.5, 0, 0}, {0, 0, 1}] *)
Ib[t_] := m / 12 * (4*0.5^2-12*0.5*(0.5+l[t])+3*1^2+12*(0.5+l[t])^2)
(* 纵摇激励力矩 *)
Mwave[t_] := L * Cos[w*t];
(* theta2 = theta1 + gamma *)
th2[t_] := th[t] + ga[t]; 
(* P点x座标 *)
xp[t_] := d*Sin[th[t]] - l[t]*Sin[th2[t]]; 
(* P点z座标 *)
zp[t_] := zg[t] - d*Cos[th[t]] + l[t]*Cos[th2[t]]; 
(* PTO系统的力 *)
Fpto[t_, c_] := -k*(l[t] - l0) - c*l'[t];
(* PTO系统的力矩 *)
Mpto[t_, crot_] := -krot*ga[t] - crot*ga'[t];

start = 200;
stop = start + 20 / w;
dim = 0.1;

P4[c_, crot_] := Sum[(c * (l'[t]^2) + crot * (ga'[t]^2)) * dim, {t, start, stop, dim}] /. NDSolve[{
    m*xp''[t] == -Fpto[t,c]*Sin[th2[t]] + Fab[t]*Cos[th2[t]], 
    m*zp''[t] == Fpto[t,c]*Cos[th2[t]] - m*g + Fab[t]*Sin[th2[t]], 
    Ib[t]*th2''[t] == Mpto[t,crot] + Fab[t]*l[t],
    M*zg''[t] == Fwave[t] + Fj[t] - M*g - Fpto[t,c]*Cos[th2[t]] - Fab[t]*Sin[th2[t]] - m0*zg''[t] - c0*zg'[t],
    Ia*th''[t] == Mwave[t] - Mpto[t,crot] + Fab[t]*d*Cos[th[t]]*Cos[th2[t]] - Fab[t]*d*Sin[th[t]]*Sin[th2[t]] - I0*th''[t] - cr0*th'[t] - Mrec*th[t],
    th'[0] == th[0] == ga[0] == ga'[0] == l'[0] == zg'[0] == 0, l[0] == l0-m/k, zg[0] == 0}, {th, ga, zg, l}, {t, start, stop}]'''

queue_task = queue.Queue()
queue_res = queue.Queue()

random.seed(1145141919)
tot = 0
threads = []
thread_num = 16
thread_running = False
thread_cmd_name = "P22"
thread_command = P22

# use_parallelize = True
use_parallelize = False

def thread_run(command, name):
    # print("thread start:", threading.current_thread())
    session = WolframLanguageSession()
    session.evaluate(command)
    # print("test exec:", session.evaluate(wlexpr("P22[90000, 0.4]")))
    while True:
        try:
            try:
                task = queue_task.get(block=True, timeout=1)
                if task is None:
                    break
                # print("got", task)
                if thread_cmd_name == "P21":
                    text = f"First[{name}[{task[0]}]]"
                else:
                    text = f"First[{name}[{task[0]}, {task[1]}]]"
                try:
                    res = session.evaluate(wlexpr(text))
                    if not isinstance(res, float):
                        raise Exception("remote error")
                    # print("res", res, threading.current_thread())
                    queue_res.put(
                        {"id": task[2], "task": task, "res": res, "ok": True})
                except Exception as e:
                    print(e, threading.current_thread())
                    queue_res.put({'id': task[2], "ok": False})
                    continue
            except queue.Empty as e:
                # print("timeout")
                if thread_running:
                    continue
                else:
                    break
        except Exception as e:
            print(e, threading.current_thread())
            queue_res.put({'id': -1, "ok": False})
            continue
    session.stop()
    if thread_running:
        print("thread exit!", threading.current_thread())


def thread_pool_init(command=P22, name="P22"):
    global tot
    global threads
    global thread_running
    global thread_cmd_name
    global thread_command
    thread_running = True
    tot = 0
    thread_cmd_name = name
    thread_command = command
    if use_parallelize:
        threads = [threading.Thread(target=thread_run, daemon=True, args=(
            command, name)) for _ in range(thread_num)]
        [t.start() for t in threads]
    else:
        session_g.evaluate(command)
    print("thread_pool_init done")

def thread_pool_exit():
    global thread_running
    global threads
    thread_running = False
    if use_parallelize:
        [t.join() for t in threads]


def parallel_evaluate(exps, retry=0, **kwargs):
    if use_parallelize:
        exps = [(*exps[i], i) for i in range(len(exps))]
        [queue_task.put(i) for i in exps]
        results = [queue_res.get(block=True) for _ in range(len(exps))]
        failed = [it['id'] for it in results if not it['ok']]
        if len(failed) > 0:
            [queue_task.put(exps[failed[i]])
            for i in range(len(failed)) if failed[i] > 0]
            results2 = [queue_res.get(block=True, timeout=15)
                        for _ in range(len(failed))]
            results = [it for it in results if it['ok']]
            results.extend(results2)
        failed = [it['id'] for it in results if not it['ok']]
        # results = [queue_res.get(block=True, timeout=15) for _ in range(len(exps))]
        results.sort(key=lambda x: x['id'])
        res = [it['res'] for it in results if it['id'] >= 0 and it['ok']]
        if len(failed) > 0 or len(results) != len(exps):
            # time.sleep(2)
            print(res, failed)
            if retry >= 3:
                raise Exception("Has failed data!!!")
            else:
                return parallel_evaluate(exps, retry=retry + 1)
        return res
    else:
        # print("exps", exps)
        if thread_cmd_name == "P21":
            task = "{i[[1]],First[%s[i[[2]]]]}" % thread_cmd_name
            values = ["{%s,%f}" % (i, exps[i][0]) for i in range(len(exps))]
        else:
            task = "{i[[1]],First[%s[i[[2]],i[[3]]]]}" % thread_cmd_name
            values = ["{%s,%f,%f}" % (i, *exps[i]) for i in range(len(exps))]
        text = "Parallelize[Table[%s, {i, {%s}}]]" % (task, ','.join(values))
        # print(text)
        results = list(session_g.evaluate(wlexpr(text)))
        # print("result", results)
        results.sort(key=lambda x: x[0])
        res = [r[1] for r in results]
        # print([isinstance(r, float) for r in res])
        all_ok = sum([0 if isinstance(r, float) else 1 for r in res]) == 0
        if not all_ok:
            print(f"retry({retry}): {exps}")
            if retry >= 3:
                print(res)
                raise Exception("Has failed data!!!")
            else:
                return parallel_evaluate(exps, retry=retry + 1)
        # print("res", res)
        return res


def func(x, y):
    global tot
    tot += 1
    return (x, y)


ERR_RAND = 0

class SA:
    def __init__(self, func, iter=thread_num, T0=thread_num, Tf=10, alpha=0.99, x_range=None, y_range=None, sx=None, sy=None, overflow=5):
        self.func = func
        self.iter = iter  # 内循环迭代次数
        self.alpha = alpha  # 降温系数
        self.T0 = T0  # 初始温度T0
        self.Tf = Tf  # 温度终值Tf
        self.T = T0  # 当前温度
        self.x_range = [0, 100000] if x_range is None else x_range
        self.y_range = [0, 1.0] if y_range is None else y_range
        self.x = [random.random() * (self.x_range[1] - self.x_range[0]) +
                  self.x_range[0] for i in range(iter)]  # 随机生成iter个x的值
        self.y = [random.random() * (self.y_range[1] - self.y_range[0]) +
                  self.y_range[0] for i in range(iter)]  # 随机生成iter个y的值
        self.thread_num = iter
        self.most_best = []
        self.history = {'f': [], 'T': []}
        self.best_F = 1e10  # 记录最好的结果
        self.best_res = None
        self.sx = sx if sx is not None else self.x[0]
        self.sy = sy if sy is not None else self.y[0]
        self.count = 0
        self.time_start = 0
        self.time_stop = 0
        self.overflow = overflow

    def generate_new(self, x, y):  # 扰动产生新解的过程
        global ERR_RAND
        while True:
            dx = self.T / self.T0 * \
                (random.random() - random.random()) * \
                (self.x_range[1] - self.x_range[0])
            dy = self.T / self.T0 * \
                (random.random() - random.random()) * \
                (self.y_range[1] - self.y_range[0])
            x_new = x + dx
            y_new = y + dy
            if (self.x_range[0] <= x_new <= self.x_range[1]) and (self.y_range[0] <= y_new <= self.y_range[1]):
                ERR_RAND = 0
                break
            else:
                # print(f"(dx, dy) = ({dx}, {dy})", f"(x_new, y_new) = ({x_new}, {y_new})",
                #     (self.x_range[0] <= x_new <= self.x_range[1]), (self.y_range[0] <= y_new <= self.y_range[1]))
                ERR_RAND += 1
                # if ERR_RAND > 30:
                #     raise Exception("ERR_RAND overflow!")
        return x_new, y_new

    def generate_directions(self):
        x, y = self.sx, self.sy

        def d(r): return self.T / self.T0 * r * (random.random() - random.random())
        dirs = [[1, 0], [0, 1], [-1, 0], [0, -1]]
        ds = [(di[0] * d(self.x_range[1] - self.x_range[0]), di[1] *
               d(self.y_range[1] - self.y_range[0])) for di in dirs]

        def apply(dd):
            # print("dd", dd, "s", (self.sx, self.sy))
            x_new, y_new = dd[0] + x, dd[1] + y
            x_new = max(x_new, self.x_range[0])
            x_new = min(x_new, self.x_range[1])
            y_new = max(y_new, self.y_range[0])
            y_new = min(y_new, self.y_range[1])
            return (x_new, y_new)
        # res_dir = [apply(dd) for dd in ds]
        res_dir = []
        dscar = 0.3
        res_dir.extend([apply((dscar*d(self.x_range[1] - self.x_range[0]), dscar *
                       d(self.y_range[1] - self.y_range[0]))) for _ in range(thread_num-len(res_dir))])
        return res_dir

    # Metropolis准则
    # 1. 如果是更好的结果则接受
    # 2. 如果不是的话则概率接受
    def Metrospolis(self, f, f_new):
        if f_new <= f:
            # print(f"accept: {f} => {f_new}")
            return True
        else:
            p = math.exp((f - f_new) / self.T)
            return random.random() < p

    def best(self):  # 获取最优目标函数值
        f_list = []  # 每次迭代之后的值
        exps = [self.func(self.x[i], self.y[i]) for i in range(self.iter)]
        results_raw = parallel_evaluate(exps, max_evaluators=self.thread_num)
        f_list = [-i for i in results_raw]
        f_best = min(f_list)
        idx = f_list.index(f_best)
        if f_best < self.best_F:
            self.best_res = (self.x[idx], self.y[idx])
        return f_best, idx  # 在该温度下，迭代L次之后目标函数的最优解和最优解的下标

    def display(self):
        plt.plot(self.history['T'], self.history['f'])
        plt.title('SA')
        plt.xlabel('T')
        plt.ylabel('f')
        plt.gca().invert_xaxis()
        plt.show()
    
    def start_run_perf(self):
        self.count = 0
        self.time_start = time.time()

    def stop_run_perf(self):
        self.time_stop = time.time()
        if self.count > 0:
            print(
                f"speed = {self.count * 60 / (self.time_stop - self.time_start)} it/min, duration = {self.time_stop - self.time_start}, count = {self.count}")

    def run_random(self):
        self.start_run_perf()
        try:
            # 外循环迭代，当前温度小于终止温度的阈值
            while self.T > self.Tf:
                # 内循环迭代
                exps = [self.func(self.x[i], self.y[i]) for i in range(self.iter)]
                exps_all = [*exps]
                xy_new = [self.generate_new(self.x[i], self.y[i])
                        for i in range(self.iter)]
                exps_new = [self.func(*xy_new[i]) for i in range(self.iter)]
                exps_all = [*exps_all, *exps_new]
                results_all = parallel_evaluate(
                    exps_all, max_evaluators=self.thread_num)
                results_all = [-i for i in results_all]
                results_new, results_raw = results_all[self.iter:], results_all[:self.iter]
                assert(len(results_raw) == len(results_new))
                for i in range(self.iter):
                    f = results_raw[i]
                    f_new = results_new[i]
                    if self.Metrospolis(f, f_new):  # 判断是否接受新值
                        self.x[i] = xy_new[i][0]  # 如果接受新值，则把新值的x,y存入x数组和y数组
                        self.y[i] = xy_new[i][1]
                # 迭代L次记录在该温度下最优解
                ft, idx = self.best()
                self.history['f'].append(ft)
                self.history['T'].append(self.T)
                # 温度按照一定的比例下降（冷却）
                self.T = self.T * self.alpha
                print("Temp now:", self.T,
                    f"F={ft}, tot={tot}, x={self.x[idx]}, y={self.y[idx]}")
                self.count += 1

            # 得到最优解
            f_best, idx = self.best()
            if f_best < self.best_F:
                print(
                    f"F={f_best}, x={self.x[idx]}, y={self.y[idx]}, count={self.count}")
            else:
                print(
                    f"F={self.best_F}, x={self.best_res[0]}, y={self.best_res[1]}, count={self.count}")
        except KeyboardInterrupt:
            thread_pool_exit()
        self.stop_run_perf()

    def run_random_climb(self):
        self.start_run_perf()
        last_f = parallel_evaluate([self.func(self.sx, self.sy)])[0]
        print("start:", f"F={last_f}, x={self.sx}, y={self.sy}")
        over = 0
        try:
            while self.T > self.Tf:
                new_pos = self.generate_directions()
                new_exps = [self.func(*p) for p in new_pos]
                new_results = parallel_evaluate(new_exps)
                # if len(new_results) != 4:
                #     raise Exception("aaaaa 4 needed!")
                self.history['f'].append(last_f)
                self.history['T'].append(self.T)
                max_index = -1
                for i in range(len(new_results)):
                    if last_f < new_results[i]:
                        max_index = i
                        last_f = new_results[i]
                        self.sx = new_pos[i][0]
                        self.sy = new_pos[i][1]
                if max_index < 0:
                    # print(f"NO WAY! last={last_f}, new_results={new_results}")
                    if over > self.overflow:
                        print("Overflows!")
                        break
                    over += 1
                    pass
                else:
                    print("update:", self.T,
                          f"F={last_f}, tot={tot}, x={self.sx}, y={self.sy}, count={self.count}")
                    over = 0
                self.T = self.T * self.alpha
                self.count += 1

            print(f"x={self.sx}, y={self.sy}, F={last_f}")
        except KeyboardInterrupt:
            thread_pool_exit()
        print("Temp now:", self.T, f"F={last_f}, tot={tot}, x={self.sx}, y={self.sy}, count={self.count}")
        self.stop_run_perf()

session_g = WolframLanguageSession() if not use_parallelize else None

if __name__ == '__main__':
    # thread_cmd_name = "P4"
    # thread_cmd_name = "P22"
    thread_cmd_name = "P21"
    if thread_cmd_name == "P4":
        # x=32241.41596050716, y=100000, F=2979.836002974773
        # Temp now: 13.7214 F=1726.2186418442523, tot=43, x=42227.56762547424, y=28799.10473477897, count=2
        # Temp now: 8.218514947124781 F=1727.2570404255514, tot=325, x=43715.00058345546, y=25654.693610991686, count=53
        if use_parallelize:
            thread_num = 6
        # sa = SA(func, x_range=[0, 100000], y_range=[0, 100000], Tf=1e-1, sx=43715.00058345546, sy=25654.693610991686)
        sa = SA(func, x_range=[0, 100000], y_range=[0, 100000], Tf=1e-1)
        thread_pool_init(command=P4, name=thread_cmd_name)
    elif thread_cmd_name == "P21":
        # x=36909.37071462201, y=0.5429676598546479, F=2295.7213258322026
        sa = SA(func, x_range=[0, 100000], y_range=[0, 1], Tf=1e-1, overflow=300)
        thread_pool_init(command=P21, name=thread_cmd_name)
    else:
        # start: F=2547.6755617509302, x=100000, y=0.44880933129114925
        # Temp now: 0.09968820536089253 F=2547.6757122465424, tot=6889, x=100000, y=0.449095895045575, count=492
        # sa = SA(func, x_range=[0, 100000], y_range=[0, 1], Tf=1e-1)
        
        # Temp now: 2.869050937515369 F=2401.167153960724, tot=2753, x=36669.47581079803, y=0, count=171
        # Temp now: 0.5631703523340716 F=2286.161260850864, tot=5345, x=36235.00330104634, y=0, count=333
        sa = SA(func, x_range=[0, 100000], y_range=[0, 1], Tf=1e-2, sx=36235.4455663941, sy=0.01)
        thread_pool_init(command=P22, name=thread_cmd_name)
    sa.run_random_climb()
    sa.display()