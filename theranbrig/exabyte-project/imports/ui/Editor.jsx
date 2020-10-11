import React, { useState } from 'react';
import styled from 'styled-components';

const Editor = ({ setElements, elements, connections, setConnections, hoveredElement }) => {
  const [xAxis, setXAxis] = useState(0);
  const [yAxis, setYAxis] = useState(0);
  const [zAxis, setZAxis] = useState(0);
  const [connection1, setConnection1] = useState(0);
  const [connection2, setConnection2] = useState(0);

  const removeFromEditor = (id) => {
    const newList = elements.filter((element) => {
      return element.id !== id;
    });
    setElements(newList);
    removeConnections(id);
  };

  createConnection = () => {
    const endPoints = elements.filter((element) => {
      console.log(element);
      console.log(connection1, connection2);
      return element.id == connection1 || element.id == connection2;
    });
    setConnections([...connections, endPoints]);
  };

  const removeConnections = (id) => {
    const filteredResults = connections.filter((connection) => {
      return connection[1].id !== id || connection[0].id === id;
    });
    setConnections(filteredResults);
  };

  return (
    <FormStyles>
      <form
        onSubmit={(e) => {
          e.preventDefault();
          setElements([
            ...elements,
            { id: elements.length + 1, x: parseInt(xAxis), y: parseInt(yAxis), z: parseInt(zAxis) },
          ]);
        }}
      >
        <label>X</label>
        <input
          type='number'
          value={xAxis}
          onChange={(e) => setXAxis(e.target.value)}
          min='-10'
          max='10'
        />
        <label>Y</label>
        <input
          type='number'
          value={yAxis}
          onChange={(e) => setYAxis(e.target.value)}
          min='-10'
          max='10'
        />
        <label>Z</label>
        <input
          type='number'
          value={zAxis}
          onChange={(e) => setZAxis(e.target.value)}
          min='-10'
          max='10'
        />
        <button type='submit'>Add</button>
      </form>
      <h3>Current Points</h3>
      <table>
        <thead>
          <tr>
            <th>no.</th>
            <th>X</th>
            <th>Y</th>
            <th>Z</th>
            <th></th>
          </tr>
        </thead>
        <tbody>
          {elements.map((element, idx) => (
            <tr key={element.id} className={element.id === hoveredElement ? 'hovered' : ''}>
              <td>{element.id}</td>
              <td>{element.x}</td> <td>{element.y}</td>
              <td>{element.z}</td>
              <td>
                <button className='remove-button' onClick={() => removeFromEditor(element.id)}>
                  X
                </button>
              </td>
            </tr>
          ))}
        </tbody>
      </table>

      <h3>Create Connection</h3>
      <p>Enter the ID of two points above to create a connection between them </p>
      <form
        onSubmit={(e) => {
          e.preventDefault();
          createConnection();
        }}
      >
        <input
          type='number'
          value={connection1}
          onChange={(e) => setConnection1(e.target.value)}
          min='1'
        />
        <input
          type='number'
          value={connection2}
          onChange={(e) => setConnection2(e.target.value)}
          min='1'
        />
        <button type='submit'>Add</button>
      </form>
    </FormStyles>
  );
};

const FormStyles = styled.div`
  width: 30%;
  display: flex;
  flex-direction: column;
  td,
  th {
    text-align: center;
    width: 20%;
  }
  tr.hovered {
    background: blue;
    color: black;
  }
  .remove-button {
    background: transparent;
    box-shadow: none;
    border: none;
    color: white;
  }
`;

export default Editor;
