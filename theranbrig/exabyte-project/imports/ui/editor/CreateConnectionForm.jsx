import React, { useState } from 'react';
import { EditorFormStyles } from './CreateNodeForm';

export const CreateConnectionForm = ({ connections, setConnections, elements }) => {
  const [connection1, setConnection1] = useState(1);
  const [connection2, setConnection2] = useState(2);

  createConnection = () => {
    const endPoints = elements.filter((element) => {
      return element.id == connection1 || element.id == connection2;
    });
    setConnections([...connections, endPoints]);
  };
  return (
    <EditorFormStyles>
      <h3>Create Connection</h3>
      <p>Enter the ID of two points above to create a connection between them </p>
      <form
        onSubmit={(e) => {
          e.preventDefault();
          createConnection();
        }}
      >
        <label>
          Start:
          <input
            type='number'
            value={connection1}
            onChange={(e) => setConnection1(e.target.value)}
            min='1'
            max={elements[elements.length - 1].id}
          />
        </label>
        <label>
          End:
          <input
            type='number'
            value={connection2}
            onChange={(e) => setConnection2(e.target.value)}
            min='1'
            max={elements[elements.length - 1].id}
          />
        </label>
        <button className='form-button' type='submit'>
          Add
        </button>
      </form>
    </EditorFormStyles>
  );
};

export default CreateConnectionForm;
