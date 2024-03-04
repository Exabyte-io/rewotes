import { VFC, useRef, useState, useEffect } from 'react';
import * as monaco from 'monaco-editor/esm/vs/editor/editor.api';
import { useLocalStorage } from '@uidotdev/usehooks';

import styles from './Editor.module.css';

// from https://github.com/microsoft/monaco-editor/blob/main/samples/browser-esm-vite-react/src/main.tsx

import Xyz from '../xyz/Xyz';
import { TwoFramesMock } from './Editor.mock';

export const Editor: VFC = () => {
    const [, shareDrawing] = useLocalStorage<Xyz>('xyzdrawing');

    const [editor, setEditor] = useState<monaco.editor.IStandaloneCodeEditor | null>(null);
    const monacoEl = useRef(null);
    const [, shareSelection] = useLocalStorage('EditorSelection', [0, 0]);

    useEffect(() => {
        if (monacoEl) {
            setEditor((editor) => {
                if (editor) return editor;

                const mEditor = monaco.editor.create(monacoEl.current!, {
                    automaticLayout: true,
                    readOnly: false,
                    value: TwoFramesMock,
                    language: 'python'
                });

                function onEditChanged() {
                    setTimeout(() => shareDrawing(Xyz.parse(mEditor.getValue())), 10);
                }

                mEditor.onDidChangeModelContent(onEditChanged);
                onEditChanged();
                mEditor.onDidChangeCursorSelection(e => {
                    shareSelection([e.selection.startLineNumber, e.selection.endLineNumber]);
                });
                // @ts-expect-error document always available in browser
                document.getElementById('import-file').addEventListener(
                    'change',
                    () => {
                        const fileInput = document.getElementById('import-file') as HTMLInputElement;
                        if (fileInput.files && fileInput.files.length > 0) {
                            // @ts-expect-error iterable in browser
                            for (const file of fileInput.files) {
                                const reader = new FileReader();
                                reader.onload = (e) => {
                                    const contents = e.target?.result as string;
                                    const text = file.name.endsWith('.poscar') ? Xyz.poscar2xyzString(contents) : contents;
                                    mEditor.setValue(mEditor.getValue() + `\n###################\n# ${file.name}\n` + text);
                                };
                                reader.readAsText(file);
                            }
                        }
                    });
                return mEditor;
            });
        }

        return () => editor?.dispose();
    }, [monacoEl.current]);

    return <div className={styles.Editor} ref={monacoEl}></div>;
};