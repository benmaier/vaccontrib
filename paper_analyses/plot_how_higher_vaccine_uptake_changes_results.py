import bfmplot as bp
pl = bp.pl
import numpy as np


N = np.array([
      9137232,
      5339517,
      46495023,
      20275029,
  ])

loc = 'de'
loc = 'en'

vacc67 = np.array([0,0.401,0.724,0.851])
vacc90 = np.array([0,0.90,0.90,0.90])

v67 = (N*vacc67).sum() / N.sum()
v90 = (N*vacc90).sum() / N.sum()

colors = [
            ['#E75740', '#F2957D'],
            ['#268D7C', '#58BDB2']
        ]

C67 = np.array([
     [0.4580029,  0.20836728],
     [0.34199717, 0.19187974]
 ])

C90 = np.array([
     [0.15779039, 0.14542971],
     [0.21461812, 0.34545187],
 ])

print(C67/C67.sum())
print(C90/C90.sum())

print(C67.sum(axis=0))
print(C90.sum(axis=0))

indices = [
            [1,1],
            [0,1],
            [1,0],
            [0,0],
          ]
labels = [ ['u → u', 'v → u'],
           ['u → v', 'v → v']
           ]

def make_bar(C,indices):

    new = []
    cols = []
    labs = []
    for n in indices:
        new.append( C[n[0],n[1]] )
        cols.append(colors[n[0]][n[1]])
        labs.append(labels[n[0]][n[1]])

    cum = np.cumsum(new)
    cum = cum[::-1]
    cols = cols[::-1]
    labs = labs[::-1]

    return cum, cols, labs


left, right = -0.65, 1.65

fig, ax = pl.subplots(1,1,figsize=(3.3,3.0))

#ax.plot([left, right], [1.2,1.2],':',c='#888888')
ax.plot([left, right], [1.,1.],':',c='#333333',zorder=-1000)
ax.text(right,1.,'R = 1',ha='right',va='bottom')
ax.plot([left, 1-0.4], [0.863,0.863],':',c='#333333',zorder=-1000)
ax.plot([left, 0-0.4], [1.200,1.200],':',c='#333333',zorder=-1000)

for i, C in enumerate([C67, C90]):
    cum, cols, lbls = make_bar(C, indices)
    for bar, color, lbl in zip(cum, cols, lbls):
        ax.bar([i],[bar],color=color)
        ax.text(i,bar-0.01,lbl,ha='center',va='top',color='w',fontsize='small')

ax.set_ylim(0,1.205)
ax.set_xlim(left,right)
ax.set_yticks([0,0.5,C90.sum(),1.0,C67.sum()])

if loc == 'de':
    ax.set_ylabel('Effektive Reproduktionszahl R')
    ax.set_xlabel('Impfquote')
else:
    ax.set_xlabel('vaccine uptake')
    ax.set_ylabel('effective reproduction number R')
ax.set_xticks([0,1])
ax.set_xticklabels([
                    '{0:2.0f}%'.format(v67*100),
                    '{0:2.0f}%'.format(v90*100),
                    ])
bp.strip_axis(ax)

fig.tight_layout()

if loc == 'de':
    fig.savefig('vaccine_uptake_de.png',dpi=150)
    fig.savefig('vaccine_uptake_de.pdf',dpi=150)
else:
    fig.savefig('vaccine_uptake_en.png',dpi=150)
    fig.savefig('vaccine_uptake_en.pdf',dpi=150)

pl.show()

