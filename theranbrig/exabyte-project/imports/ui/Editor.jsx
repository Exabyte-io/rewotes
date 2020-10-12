import React, { useState } from 'react';
import styled from 'styled-components';
import CurrentPointsTable from './CurrentPointsTable';
import AddCoordinateForm from './CreateNodeForm';
import CreateConnectionForm from './CreateConnectionForm';
import SaveMaterial from './SaveMaterial';

const Editor = ({ setElements, elements, connections, setConnections, hoveredElement }) => {
  const [name, setName] = useState('');

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
      <h1>Materials Editor</h1>
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
          <SaveMaterial elements={elements} connections={setConnections} />
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
