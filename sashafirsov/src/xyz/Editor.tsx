import { VFC, useRef, useState, useEffect } from 'react';
import * as monaco from 'monaco-editor/esm/vs/editor/editor.api';
import styles from './Editor.module.css';

// from https://github.com/microsoft/monaco-editor/blob/main/samples/browser-esm-vite-react/src/main.tsx

export const Editor: VFC = () => {
    const [editor, setEditor] = useState<monaco.editor.IStandaloneCodeEditor | null>(null);
    const monacoEl = useRef(null);

    useEffect(() => {
        if (monacoEl) {
            setEditor((editor) => {
                if (editor) return editor;

                return monaco.editor.create(monacoEl.current!, {
                    automaticLayout: true,
                    value: ['function x() {', '\tconsole.log("Hello world!");', '}'].join('\n'),
                    language: 'typescript'
                });
            });
        }

        return () => editor?.dispose();
    }, [monacoEl.current]);

    return <div className={styles.Editor} ref={monacoEl}></div>;
};