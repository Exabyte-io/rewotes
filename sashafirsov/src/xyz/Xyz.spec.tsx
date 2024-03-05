import { render } from '@testing-library/react';

import Xyz from './Xyz';
import { poscarMgO, TwoFramesMock } from './Editor.mock';

describe('Xyz', () => {
    it('parse', () => {
        const xyz = Xyz.parse(TwoFramesMock)
        expect(xyz.slides.length).toEqual(2);
        expect(xyz.slides[0].comment.trim()).toEqual('Created by chemcoord http://chemcoord.readthedocs.io/en/latest/');
        expect(xyz.slides[1].comment.trim()).toEqual('# pyridine molecule');

        expect(xyz.slides[0].elements.length).toEqual(6);
        expect(xyz.slides[1].elements.length).toEqual(11);

        expect(xyz.slides[0].elements[5]).deep.equal({
            "element": "H",
            "x": 3.260455,
            "y": 0.5,
            "z": -0.872893,
            "sourceLine": 8 // -1 as it is 0-based
        });
        expect(xyz.slides[1].elements[10]).deep.equal({
            "element": "H",
            "x": -0.180226841,
            "y": -1.796059882,
            "z": -0.917077970,
            "sourceLine": 21 // -1 as it is 0-based
        });
    });
    it('poscar2xyzString', () => {
        const xyzStr = Xyz.poscar2xyzString(poscarMgO)
        expect(xyzStr).toEqual(
`# scaling 1.0
2
MgO Fm-3m (No. 225)
Mg  0.000000 0.000000 0.000000 Mg
O  0.500000 0.500000 0.500000 O`);
    });
});
