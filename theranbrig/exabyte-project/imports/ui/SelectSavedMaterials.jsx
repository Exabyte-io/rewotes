import React from 'react';
import styled from 'styled-components';

const SelectSavedMaterials = ({ selectedMaterial, selectMaterial, materials }) => {
  const onChangeHandler = (id) => {
    selectMaterial(id);
  };
  return (
    <SelectStyles>
      <h3>Select Saved Material</h3>
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
    </SelectStyles>
  );
};

export const SelectStyles = styled.div`
  select {
    padding: 20px;
    background: #273746;
    border-radius: 10px;
    border: none;
    color: white;
    box-shadow: none;
    width: 100%;
    margin-top: 20px;
  }
`;

export default SelectSavedMaterials;
