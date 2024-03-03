export type XyzArgs = [string, number, number, number];

export class ElementXyz {
    element: string;
    x: number;
    y: number;
    z: number;

    constructor(...args: XyzArgs | [string]) {
        if (args.length === 1 && typeof args[0] === 'string') {
            const [element, x, y, z] = args[0].split(/\s+/);
            this.element = element;
            this.x = parseFloat(x);
            this.y = parseFloat(y);
            this.z = parseFloat(z);
        } else {
            const [element, x, y, z] = args as XyzArgs;
            this.element = element;
            this.x = x;
            this.y = y;
            this.z = z;
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