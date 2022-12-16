import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import csv

data = open('Implicit_Localization_out.txt', 'r').read()
cnt = 0
mas = []
tmp_str = ''

for i in data:
    if i == '\n':
        cnt += 1
    tmp_str = tmp_str + i

    if cnt == 252:
        cnt = 0
        mas.append(tmp_str[:-1])
        tmp_str = ''
print(mas)

data2 = open('Localization_exact.txt', 'r').read()
cnt = 0
mas2 = []
tmp_str = ''

for i in data:
    if i == '\n':
        cnt += 1
    tmp_str = tmp_str + i

    if cnt == 252:
        cnt = 0
        mas2.append(tmp_str[:-1])
        tmp_str = ''
print(mas2)

X1 = []
Y1 = []
with open('Localization_start.txt', 'r') as datafile:
    plotting = csv.reader(datafile, delimiter=' ')

    for ROWS in plotting:
        X1.append(float(ROWS[0]))
        Y1.append(float(ROWS[1]))

fig, ax1 = plt.subplots(1, 1)
plt.ylim([0,30])
plt.xlim([-15,15])

anim_cnt = 0
def animate(i):
    global anim_cnt
    try:
        tmp_data = mas[anim_cnt]
        tmp_data2 = mas2[anim_cnt]
        if tmp_data != '' and tmp_data2 != '':
            tmp_lines = tmp_data.split('\n')
            tmp_lines2 = tmp_data2.split('\n')
            xs = []
            ys = []
            xs2 = []
            ys2 = []
            for lines in tmp_lines:
                x, y = lines.split(' ')
                xs.append(float(x))
                ys.append(float(y))

            for lines in tmp_lines2:
                x, y = lines.split(' ')
                xs2.append(float(x))
                ys2.append(float(y)) 
            ax1.clear()
            
            ax1.plot(xs2, ys2, color = 'b', label='exact')
            ax1.plot(xs, ys, color = 'g', label='solution')
            ax1.plot(X1, Y1, color = 'r', label='start')
            ax1.legend()
            
            anim_cnt += 1
    except Exception as e:
        time.sleep(1000)

ani = animation.FuncAnimation(fig, animate, interval=200)
plt.show()


