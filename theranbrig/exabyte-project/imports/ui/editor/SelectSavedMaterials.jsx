import React from 'react';
import styled from 'styled-components';

const SelectSavedMaterials = ({ selectedMaterial, selectMaterial, materials }) => {
  const onChangeHandler = (id) => {
    selectMaterial(id);
  };
  return (
    <SelectStyles>
      <div className='select-div'>
        <label>Select Saved Material</label>
        <div className='select-container'>
          <select value={selectedMaterial} onChange={(e) => onChangeHandler(e.target.value)}>
            <option value={null}>-</option>
            {materials
              .filter((material) => material.connections && material.name && material.nodes)
              .map((material) => (
                <option key={material._id} value={material._id}>
                  {material.name}
                </option>
              ))}
          </select>
        </div>
      </div>
    </SelectStyles>
  );
};

export const SelectStyles = styled.div`
  .select-div {
    position: relative;
    width: 100%;
  }
  select::-ms-expand {
    display: none;
  }
  .select-container {
    background: #273746;
    padding: 20px;
    height: 29px;
    margin-top: 10px;
    border-radius: 10px;
  }
  .select-div:after {
    content: '<>';
    font: 17px 'Consolas', monospace;
    color: white;
    -webkit-transform: rotate(90deg);
    -moz-transform: rotate(90deg);
    -ms-transform: rotate(90deg);
    transform: rotate(90deg);
    right: 30px;
    top: 50px;
    padding: 0 0 2px;
    border-bottom: 1px solid white;
    position: absolute;
    pointer-events: none;
  }

  .select-div select {
    -webkit-appearance: none;
    -moz-appearance: none;
    appearance: none;
    width: 100%;
    float: right;
    font-size: 16px;
    line-height: 1.75;
    color: white;
    padding: 0 10px;
    background-color: transparent;
    background-image: none;
    border: none;
    border-bottom: 1px solid white;
    -ms-word-break: normal;
    word-break: normal;
  }
`;

export default SelectSavedMaterials;
