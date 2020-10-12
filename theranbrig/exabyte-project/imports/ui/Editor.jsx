import React from 'react';
import styled from 'styled-components';
import CurrentPointsTable from './CurrentPointsTable';
import AddCoordinateForm from './CreateNodeForm';
import CreateConnectionForm from './CreateConnectionForm';
import SaveMaterial from './SaveMaterial';
import SelectSavedMaterials from './SelectSavedMaterials';

const Editor = ({
  setElements,
  elements,
  connections,
  setConnections,
  hoveredElement,
  materials,
  selectMaterial,
  name,
  setName,
  selectedMaterial,
  clearEditor,
}) => {
  const removeFromEditor = (id) => {
    if (window.confirm('Are you sure you want to remove this?')) {
      const newList = elements.filter((element) => {
        return element.id !== id;
      });
      setElements(newList);
      removeConnections(id);
    }
  };

  const removeConnections = (id) => {
    const filteredResults = connections.filter((connection) => {
      return connection[1].id !== id && connection[0].id !== id;
    });
    setConnections(filteredResults);
  };

  return (
    <FormStyles>
      <h1>{selectedMaterial ? `${selectedMaterial.name} Editor` : 'Materials Editor'}</h1>
      <SelectSavedMaterials materials={materials} selectMaterial={selectMaterial} />
      <AddCoordinateForm elements={elements} setElements={setElements} />
      {elements.length ? (
        <>
          <CurrentPointsTable
            elements={elements}
            hoveredElement={hoveredElement}
            removeFromEditor={removeFromEditor}
          />
          {elements.length > 1 ? (
            <CreateConnectionForm
              connections={connections}
              setConnections={setConnections}
              elements={elements}
            />
          ) : null}
          <SaveMaterial
            elements={elements}
            connections={connections}
            name={name}
            setName={setName}
            clearEditor={clearEditor}
            selectedMaterial={selectedMaterial}
          />
        </>
      ) : null}
    </FormStyles>
  );
};

const FormStyles = styled.div`
  width: 30%;
  display: flex;
  flex-direction: column;
  padding: 10px;
  overflow-y: scroll;
  h3,
  h1 {
    font-weight: 300;
    margin: 30px 0 5px;
  }
`;

export default Editor;
