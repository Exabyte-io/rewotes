import React, { useState } from 'react';
import styled from 'styled-components';

const CreateNodeForm = ({ elements, setElements }) => {
  const [xAxis, setXAxis] = useState(0);
  const [yAxis, setYAxis] = useState(0);
  const [zAxis, setZAxis] = useState(0);
  return (
    <EditorFormStyles>
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
    </EditorFormStyles>
  );
};

export const EditorFormStyles = styled.div`
  form {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    background: #273746;
    padding: 20px;
    border-radius: 10px;
    margin-top: 20px;
  }

  .form-button {
    padding: 5px 15px;
    border-radius: 9999rem;
    box-shadow: none;
    border: none;
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
  .remove-button {
    background: transparent;
    border: 1px solid white;
    padding: 5px 15px;
    box-shadow: none;
    border-radius: 9999rem;
    margin-top: 10px;
    img {
      height: 20px;
    }
  }
`;
export default CreateNodeForm;
