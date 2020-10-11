import React, { useState } from 'react';
import CanvasArea from './Canvas';
import Editor from './Editor';
import styled from 'styled-components';

export const App = () => {
  const [elements, setElements] = useState([]);
  const [connections, setConnections] = useState([]);
  const [hoveredElement, setHoveredElement] = useState(null);
  return (
    <MainContentStyling>
      <Editor
        elements={elements}
        setElements={setElements}
        connections={connections}
        setConnections={setConnections}
        hoveredElement={hoveredElement}
      />
      <CanvasArea
        elements={elements}
        connections={connections}
        setHoveredElement={setHoveredElement}
      />
    </MainContentStyling>
  );
};

const MainContentStyling = styled.div`
  background: #1c2833;
  color: white;
  display: flex;
  flex-direction: row;
  height: 100vh;
  margin: 0;
  padding: 0;
  font-weight: 300;
  overflow: hidden;
`;
