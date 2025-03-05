import React, { useState } from 'react';
import CanvasArea from './canvas/Canvas';
import Editor from './editor/Editor';
import styled from 'styled-components';
import { useTracker } from 'meteor/react-meteor-data';
import { Materials } from '/imports/api/materials';
export const App = () => {
  const [elements, setElements] = useState([]);
  const [connections, setConnections] = useState([]);
  const [hoveredElement, setHoveredElement] = useState(null);
  const [selectedMaterial, setSelectedMaterial] = useState(null);
  const [name, setName] = useState('');

  const { materials } = useTracker(() => {
    Meteor.subscribe('materials');
    return { materials: Materials.find().fetch() };
  });

  const clearEditor = (selected) => {
    setSelectedMaterial(selected ? selected : null);
    setConnections(selected ? [...selected.connections] : []);
    setElements(selected ? [...selected.nodes] : []);
    setName(selected ? selected.name : '');
  };

  const selectMaterial = (id) => {
    const selected = materials.filter((material) => material._id === id)[0];
    clearEditor(selected);
  };

  return (
    <MainContentStyling>
      <Editor
        elements={elements}
        setElements={setElements}
        connections={connections}
        setConnections={setConnections}
        hoveredElement={hoveredElement}
        materials={materials}
        selectMaterial={selectMaterial}
        name={name}
        setName={setName}
        selectedMaterial={selectedMaterial}
        clearEditor={clearEditor}
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
