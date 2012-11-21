from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     AnswerTestingTest, \
     sim_dir_load
from yt.frontends.enzo.answer_testing_support import \
     requires_outputlog

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)
_fields = ('Density', 'Pressure')

class TestFryxellImage(AnswerTestingTest):
    _type_name = "fryxell_image"
    _attrs = ("field", )

    def __init__(self, pf, field):
        self.pf = pf
        self.field = field

    def run(self):
        sl = self.pf.h.slice(2, 0.5)
        frb = FixedResolutionBuffer(sl, (0.0, 1.0, 0.0, 1.0), (200,200))
        dd = frb[self.field]
        return np.array([dd.mean(), dd.std(), dd.min(), dd.max()])

    def compare(self, new_result, old_result):
        assert_allclose(new_result, old_result, rtol=1e-13, atol=0)

@requires_outputlog(_dir_name, _pf_name)
def test_fryxell():
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)
    sim.get_time_series()
    for pf in sim:
        for field in _fields:
            yield TestFryxellImage(pf, field)
