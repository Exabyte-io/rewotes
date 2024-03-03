export type XyzArgs = [string, number, number, number,number];

export class ElementXyz {
    element: string;
    x: number;
    y: number;
    z: number;
    sourceLine:number;

    constructor(...args: XyzArgs | [string]) {
        if (args.length === 1 && typeof args[0] === 'string') {
            const [element, x, y, z,sourceLine] = args[0].split(/\s+/);
            this.element = element;
            this.x = parseFloat(x);
            this.y = parseFloat(y);
            this.z = parseFloat(z);
            this.sourceLine = parseInt(sourceLine)
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
}