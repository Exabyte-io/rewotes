import { Card, ControlGroup, TextArea, InputGroup, Tag, Button } from '@blueprintjs/core'
import React, { useEffect } from 'react'
import styles from './SourceEditor.module.scss';
import { useSourceContext } from '../../context/SourceContext';
import { useAtomsContext, Atom } from '../../context/AtomsContext';
import VStack from '../../components/utils/VStack';
import ViewHeading from '../../components/view_heading/ViewHeading';
import { Vector3 } from 'three';

function SourceEditor() {
  const { source, setSource, importSource, sourceName, setSourceName, isValidXYZFormat, setIsValidXYZFormat, saveSourceToLocalStorage } = useSourceContext();
  const { atoms, updateAtoms } = useAtomsContext();

  type ParsedXYZResult = {
    isValid: boolean;
    atoms?: Atom[];
  };

  const parseXYZ = (text: string): ParsedXYZResult => {
    const FIRST_TWO_LINES = 2;
    const lines = text.split('\n');
    const atomCount = parseInt(lines[0].trim(), 10);
    const newAtoms: Atom[] = [];

    if (isNaN(atomCount) || lines.length < atomCount + FIRST_TWO_LINES) {
      return { isValid: false }; //Must be the exact amount of lines 
    }

    for (let i = FIRST_TWO_LINES; i < atomCount + FIRST_TWO_LINES; i++) {
      const line = lines[i].split(/\s+/);
      if (line.length !== 4) { // "element, x, y, z" -- exactly 4 objects
        return { isValid: false };
      }
      const element = line[0]; //TODO: Add validation against elements enum

      for (let j = 1; j < line.length; j++) {
        const elementValue = line[j];
        if (!/^-?\d*\.?\d+$/.test(elementValue)) {
          return { isValid: false };
        }
      }

      const x = parseFloat(line[1]);
      const y = parseFloat(line[2]);
      const z = parseFloat(line[3]);
      if (isNaN(x) || isNaN(y) || isNaN(z)) {
        return { isValid: false };
      }

      newAtoms.push({
        id: i - FIRST_TWO_LINES,
        element: element,
        position: new Vector3(x, y, z),
      });
    }

    return { isValid: true, atoms: newAtoms };
  };

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
    else console.log('problem at handle chnge');
  };

  const handleSave = () => {
    //TODO: save to local storage 
  }

  return (
    <Card className={styles.SourceEditor}>
      <VStack>
        <ViewHeading>
          <h4>Source editor</h4>
        </ViewHeading>

        <ControlGroup>
          <InputGroup readOnly={true} value={sourceName} fill={true} rightElement={<Button text='Save' onClick={handleSave} />} />
          <Tag
            style={{ width: '30%' }}
            intent={isValidXYZFormat ? 'success' : 'warning'}
          >
            {isValidXYZFormat ? 'Correct XYZ pattern' : 'Wrong XYZ pattern'}
          </Tag>
        </ControlGroup>
        <TextArea
          id='text-editor'
          className={styles.TextEditor}
          value={source}
          onChange={handleChange}
        />
      </VStack>
    </Card>
  )
}

export default SourceEditor

