import React, { useState } from 'react';
import { EditorFormStyles } from './CreateNodeForm';
import Materials from '/imports/api/materials';

const SaveMaterial = ({ elements, connections, name, setName }) => {
  const handleSubmit = (e) => {
    e.preventDefault();
    if (!name) return;
    Meteor.call('materials.insert', name, elements, connections);
  };

  return (
    <EditorFormStyles>
      <h3>Save Current Model</h3>
      <form onSubmit={handleSubmit}>
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
    </EditorFormStyles>
  );
};

export default SaveMaterial;
