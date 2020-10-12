import React from 'react';
import { EditorFormStyles } from './CreateNodeForm';

 const SaveMaterial = () => {
  return (
    <EditorFormStyles>
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
    </EditorFormStyles>
  );
};

export default SaveMaterial;
