# This file is part of TetWild, a software for generating tetrahedral meshes.
#
# Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu> and the TetWild contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
#
# Created by trelau on 07/06/2018.
import os
import subprocess
import time
import unittest
from shutil import copyfile


class TestTetWild(unittest.TestCase):
    """
    Test cases for TetWild.
    """

    @classmethod
    def setUpClass(cls):
        if not os.path.exists('./scratch'):
            os.mkdir('scratch')

    @classmethod
    def run_tetwild(cls, fin, fout):
        """
        Run the TetWild executable on the input file.

        :param str fin: The input file.
        :param str fout: The output file name.

        :return: The runtime.
        :rtype: float
        """
        copyfile('./models/{}'.format(fin), './scratch/{}'.format(fin))
        arg = '--input=./scratch/{} --output=./scratch/{}'.format(fin, fout)
        cmd = ' '.join(['TetWild', arg])
        start = time.time()
        subprocess.call(cmd)
        return time.time() - start

    def test_12246(self):
        """
        Thingi10k: https://ten-thousand-models.appspot.com/detail.html?file_id=39549
        """
        runtime = self.run_tetwild('12246.stl', '12246.msh')
        self.assertLessEqual(runtime, 130.)
