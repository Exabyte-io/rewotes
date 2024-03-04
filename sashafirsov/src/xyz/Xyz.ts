export type XyzArgs = [string, number, number, number, number];

export class ElementXyz {
    element: string;
    x: number;
    y: number;
    z: number;
    sourceLine: number;

    constructor(...args: XyzArgs | [string]) {
        if (args.length === 1 && typeof args[0] === 'string') {
            const [element, x, y, z, sourceLine] = args[0].split(/\s+/);
            this.element = element;
            this.x = parseFloat(x);
            this.y = parseFloat(y);
            this.z = parseFloat(z);
            this.sourceLine = parseInt(sourceLine);
        } else {
            const [element, x, y, z, sourceLine] = args as XyzArgs;
            this.element = element;
            this.x = x;
            this.y = y;
            this.z = z;
            this.sourceLine = sourceLine;
        }
    }
}

export class XyzSlide {
    comment = '';
    elements: ElementXyz[] = [];
}

export default class Xyz {
    name = '';
    slides: XyzSlide[] = [];

    addSlide(xyzSlide: XyzSlide) {
        this.slides.push(xyzSlide);
    }

    static parse(xyzString: string) {
        const xyz = new Xyz();

        const lines = xyzString.split('\n');

        for (let i = 0; i < lines.length; i++) {
            const line = lines[i].trim();
            if (!line || line.startsWith('#')) // ignore blank lines and comments
                continue;

            const n = parseInt(lines[i++]);
            const comment = lines[i++];
            const elements = lines.slice(i, i + n).map((line, k) => {
                const [element, x, y, z] = line.split(/\s+/);
                return new ElementXyz(element, parseFloat(x), parseFloat(y), parseFloat(z), i + k + 1);
            });
            const slide = new XyzSlide();
            slide.comment = comment;
            slide.elements = elements;
            xyz.addSlide(slide);
            i += n - 1;
        }
        return xyz;
    }

    static poscar2xyzString(contents: string) {
        const [comment, scaling, _a1, _a2, _a3, ionSpecies, ionNumbers, _directOrCartesian, ...positions] = contents.split('\n');

        const species = ionSpecies.trim().split(/\s+/);
        const speciesCount = ionNumbers.trim().split(/\s+/).map(n => parseInt(n));
        const xyzLines = [
            '# scaling ' + scaling,
            speciesCount.reduce((total, v) => total + v, 0),
            comment
        ];
        let pi = 0;
        for (let si = 0; si < species.length; si++) {
            for (let sci = 0; sci < speciesCount[si]; sci++, pi++)
                xyzLines.push(`${species[si]} ${positions[pi]}`);
        }
        return xyzLines.join('\n');
    }
}