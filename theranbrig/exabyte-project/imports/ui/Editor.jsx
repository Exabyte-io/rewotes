import React, { useState } from 'react';
import styled from 'styled-components';

const Editor = ({ setElements, elements, connections, setConnections, hoveredElement }) => {
  const [xAxis, setXAxis] = useState(0);
  const [yAxis, setYAxis] = useState(0);
  const [zAxis, setZAxis] = useState(0);
  const [connection1, setConnection1] = useState(1);
  const [connection2, setConnection2] = useState(2);
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

  createConnection = () => {
    const endPoints = elements.filter((element) => {
      console.log(element);
      console.log(connection1, connection2);
      return element.id == connection1 || element.id == connection2;
    });
    setConnections([...connections, endPoints]);
  };

  const removeConnections = (id) => {
    const filteredResults = connections
      .filter((connection) => {
        return connection[1].id !== id;
      })
      .filter((connection) => {
        return connection[0].id !== id;
      });
    setConnections(filteredResults);
  };

  return (
    <FormStyles>
      <h1>Materials Editor</h1>
      <h3>Add Point to Grid</h3>
      <p>Enter the coordinates below</p>
      <form
        className='axis-form'
        onSubmit={(e) => {
          e.preventDefault();
          setElements([
            ...elements,
            {
              id: !elements.length ? 1 : elements[elements.length - 1].id + 1,
              x: parseInt(xAxis),
              y: parseInt(yAxis),
              z: parseInt(zAxis),
            },
          ]);
        }}
      >
        <label>
          X:
          <input
            type='number'
            value={xAxis}
            onChange={(e) => setXAxis(e.target.value)}
            min='-10'
            max='10'
          />
        </label>
        <label>
          Y:
          <input
            type='number'
            value={yAxis}
            onChange={(e) => setYAxis(e.target.value)}
            min='0'
            max='10'
          />
        </label>
        <label>
          Z:
          <input
            type='number'
            value={zAxis}
            onChange={(e) => setZAxis(e.target.value)}
            min='-10'
            max='10'
          />
        </label>
        <button className='form-button' type='submit'>
          Add
        </button>
      </form>
      <h3>Current Points</h3>
      <div>
        <table className='table-container'>
          <thead>
            <tr>
              <th>ID</th>
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
      </div>
      {elements.length > 1 ? (
        <>
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
        </>
      ) : null}
      <h3>Save Current Model</h3>
      <form>
        <label>
          Material Name:
          <input
            style={{ textAlign: 'left' }}
            type='text'
            value={name}
            onChange={(e) => setName(e.target.value)}
          />
        </label>
        <button className='form-button' type='submit'>
          Save
        </button>
      </form>
    </FormStyles>
  );
};

const FormStyles = styled.div`
  width: 30%;
  display: flex;
  flex-direction: column;
  padding: 10px;
  overflow-y: scroll;
  td,
  th {
    text-align: center;
    width: 20%;
    padding: 10px 0;
  }
  table {
    margin-top: 20px;
    border-spacing: 0px;
    width: 100%;
  }
  thead {
    display: table-header-group;
    vertical-align: middle;
  }
  th {
    background: #273746;
    font-weight: 300;
  }
  tr:nth-child(even) {
    background: #808b96;
  }
  tr:nth-child(odd) {
    background: #566573;
  }

  tr.hovered {
    background: #2471a3;
    color: white;
  }
  .table-container {
    width: 100%;
    padding: 10px;
    border-radius: 10px;
    background: #273746;
  }
  .remove-button {
    background: transparent;
    box-shadow: none;
    border: none;
    color: white;
  }
  input {
    background: transparent;
    border: none;
    border-bottom: 1px solid white;
    box-shadow: none;
    margin-left: 5px;
    color: white;
    text-align: center;
  }
  form {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    background: #273746;
    padding: 20px;
    border-radius: 10px;
  }

  .form-button {
    padding: 5px 15px;
    border-radius: 9999rem;
    box-shadow: none;
    border: none;
  }
  h3,
  h1 {
    font-weight: 300;
    margin: 30px 0 5px;
  }
`;

export default Editor;
