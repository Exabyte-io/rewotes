import { Card, ControlGroup, TextArea, InputGroup, Tag, Button } from '@blueprintjs/core'
import React, { useEffect } from 'react'
import styles from './SourceEditor.module.scss';
import { useSourceContext } from '../../context/SourceContext';
import { useAtomsContext, Atom } from '../../context/AtomsContext';
import VStack from '../../components/utils/VStack';
import ViewHeading from '../../components/view_heading/ViewHeading';
import { Vector3 } from 'three';
import { poscarToXYZ } from '../../actions/poscarToXyz';
import parseXYZ from '../../actions/parseXYZ';

function SourceEditor() {
  const { source, setSource, sourceName, isValidXYZFormat, setIsValidXYZFormat } = useSourceContext();
  const { updateAtoms } = useAtomsContext();

  useEffect(() => {
    const parsedResult = parseXYZ(source);
    setIsValidXYZFormat(parsedResult.isValid);
  }
  );

  const handleChange = (event: React.ChangeEvent<HTMLTextAreaElement>) => {
    const newValue = event.target.value;
    setSource(newValue);
    const parsedResult = parseXYZ(newValue);
    setIsValidXYZFormat(parsedResult.isValid);
    if (parsedResult.isValid && parsedResult.atoms) {
      parsedResult.atoms.forEach((atom) => {
        updateAtoms(atom.id, atom.position);
        console.log('updating atoms');
      });
    }
    else console.log('problem at handle change');
  };

  const handleConvert = () => {
    if (!isValidXYZFormat) {
      const data = poscarToXYZ(source);
      if (data) setSource(data); 
    }
    else throw new Error('Failed to convert: already a XYZ format');
  }

  return (
    <Card className={styles.SourceEditor}>
      <VStack>
        <ViewHeading>
          <h4>Source editor</h4>
        </ViewHeading>

        <ControlGroup data-testid='control-group'>
          <InputGroup readOnly={true} value={sourceName} fill={true} rightElement={<Button text='Convert' onClick={handleConvert} />} />
          <Tag
            data-testid='xyz-validity-tag'
            style={{ width: '30%', textAlign: 'center' }}
            intent={isValidXYZFormat ? 'success' : 'warning'}
          >
            {isValidXYZFormat ? 'Correct' : 'Wrong'} XYZ format
          </Tag>
        </ControlGroup>
        <TextArea
          id='text-editor'
          aria-label='source-editor-textarea'
          className={styles.TextEditor}
          value={source}
          onChange={handleChange}
        />
      </VStack>
    </Card>
  )
}

export default SourceEditor

