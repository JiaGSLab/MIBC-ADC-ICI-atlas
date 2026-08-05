from pathlib import Path
import importlib.util, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text

ROOT=Path('/data/home/chenwenhui/fig6')
OUT=ROOT/'figures_vector/main/candidate_6f_corrected_axes_width80_20260728'
FONT=ROOT/'figures_vector/main/candidate_compact_20260728/fonts'
SRC=Path(__file__).with_name('plot_fig6_revised_panels.py')
spec=importlib.util.spec_from_file_location('rf', SRC); rf=importlib.util.module_from_spec(spec); sys.modules['rf']=rf; spec.loader.exec_module(rf)
for p in FONT.glob('*.ttf'): font_manager.fontManager.addfont(str(p))
plt.rcParams.update({'font.family':'Arial','font.sans-serif':['Arial'],'pdf.fonttype':42,'ps.fonttype':42})
ladder,diffs=rf.load_composite_ladder()
expected={'celltype_proportion':0.439628,'CN_proportion':0.715170,'CN_coupling':0.845201}
for r in ladder.itertuples():
    assert abs(r.auc-expected[r.model_id]) < 1e-6, (r.model_id,r.auc)
OUT.mkdir(parents=True,exist_ok=True)
fig,ax=plt.subplots(figsize=(3.780,1.960))
rf.draw_ladder(ax,ladder,diffs)
# Width only: native uncompressed plot box was 0.59 of figure width; retain height
# and reduce its physical width to 80% (0.472), keeping coordinates unchanged.
ax.set_position([.50,.30,.472,.54])
for t in fig.findobj(Text):
    if t.get_visible(): t.set_fontfamily('Arial'); t.set_fontsize(8)
fig.savefig(OUT/'fig6f_auc_ladder_corrected_axes_width80.pdf',format='pdf',facecolor='white')
