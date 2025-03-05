import React from 'react';
import SplitPane from 'split-pane-react';
import useResizeable from '../../hooks/useResizeable';

const ResizablePane = ({ children }) => {
    // Style states
    const { sizes, handleSizeChange } = useResizeable();

    return (
        <SplitPane split='vertical' sizes={sizes} onChange={handleSizeChange}>
            {children}
        </SplitPane>
    );
};

export default ResizablePane;
