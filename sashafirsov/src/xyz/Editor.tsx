import { VFC, useRef, useState, useEffect } from 'react';
import * as monaco from 'monaco-editor/esm/vs/editor/editor.api';
import { useLocalStorage } from '@uidotdev/usehooks';

import styles from './Editor.module.css';

// from https://github.com/microsoft/monaco-editor/blob/main/samples/browser-esm-vite-react/src/main.tsx

import Xyz, { ElementXyz, XyzSlide } from '../xyz/Xyz';

export const Editor: VFC = () => {
    const [drawing, saveDrawing] = useLocalStorage<Xyz>('xyzdrawing');

    const [editor, setEditor] = useState<monaco.editor.IStandaloneCodeEditor | null>(null);
    const monacoEl = useRef(null);

    useEffect(() => {
        if (monacoEl) {
            setEditor((editor) => {
                if (editor) return editor;

                const mEditor = monaco.editor.create(monacoEl.current!, {
                    automaticLayout: true,
                    readOnly: false,
                    value:
                        `6
Created by chemcoord http://chemcoord.readthedocs.io/en/latest/
O 0.000000 0.000000  0.000000
H 0.758602 0.000000  0.504284
H 0.260455 0.000000 -0.872893
O 3.000000 0.500000  0.000000
H 3.758602 0.500000  0.504284
H 3.260455 0.500000 -0.872893
11
# pyridine molecule 
C       -0.180226841      0.360945118     -1.120304970
C       -0.180226841      1.559292118     -0.407860970
C       -0.180226841      1.503191118      0.986935030
N       -0.180226841      0.360945118      1.29018350
C       -0.180226841     -0.781300882      0.986935030
C       -0.180226841     -0.837401882     -0.407860970
H       -0.180226841      0.360945118     -2.206546970
H       -0.180226841      2.517950118     -0.917077970
H       -0.180226841      2.421289118      1.572099030
H       -0.180226841     -1.699398882      1.572099030
H       -0.180226841     -1.796059882     -0.917077970
`,
                    language: 'python'
                });
                mEditor.onDidChangeModelContent(function(e) {
                    // const m = mEditor.getModel();
                    const xyz = new Xyz();

                    const lines = mEditor.getValue().split('\n').filter(t => t.trim());

                    for (let i = 0; i < lines.length; i++) {
                        const n = parseInt(lines[i++]);
                        const comment = lines[i++];
                        const elements = lines.slice(i, i + n).map((line) => {
                            const [element, x, y, z] = line.split(/\s+/);
                            return new ElementXyz(element, parseFloat(x), parseFloat(y), parseFloat(z));
                        });
                        const slide = new XyzSlide();
                        slide.comment = comment;
                        slide.elements = elements;
                        xyz.addSlide(slide);
                        i += n;
                    }
                    saveDrawing(xyz);

                });
                return mEditor;
            });
        }

        return () => editor?.dispose();
    }, [monacoEl.current]);

    return <div className={styles.Editor} ref={monacoEl}></div>;
};