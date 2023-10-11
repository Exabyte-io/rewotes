import React, { useState } from 'react';
import { EditorFormStyles } from './CreateNodeForm';
import Materials from '/imports/api/materials';

const SaveMaterial = ({ elements, connections, name, setName, selectedMaterial, clearEditor }) => {
  const handleSave = (e) => {
    e.preventDefault();
    if (!name) return;
    if (!selectedMaterial) {
      Meteor.call('materials.insert', name, elements, connections);
    } else {
      Meteor.call('materials.update', selectedMaterial._id, name, elements, connections);
    }
  };

  const handleRemove = () => {
    console.log(selectedMaterial);
    if (window.confirm('Are you sure you want to delete this material?')) {
      Meteor.call('materials.remove', selectedMaterial._id);
      clearEditor(null);
    }
  };

  return (
    <EditorFormStyles>
      <h3>Save Current Model</h3>
      <form onSubmit={handleSave}>
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
      <button onClick={handleRemove} className='remove-button' aria-label='delete material'>
        <img src='/images/delete.svg' alt='trash can' />
      </button>
    </EditorFormStyles>
  );
};

export default SaveMaterial;
