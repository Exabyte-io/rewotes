import { Card, ControlGroup, TextArea, InputGroup, Tag } from '@blueprintjs/core'
import React, { useState } from 'react'
import styles from './SourceEditor.module.scss';
import { useSourceContext } from '../../context/SourceContext';
import VStack from '../../components/utils/VStack';

function SourceEditor() {
  const { source, setSource, importSource, sourceName, setSourceName, isValidXYZFormat, setIsValidXYZFormat } = useSourceContext();

  const isValidXYZ = (content: string) => {
    const regexXYZ = /^\d+\n\s*\n(?:\w+\s+[\d.-]+\s+[\d.-]+\s+[\d.-]+\n)+$/;
    return regexXYZ.test(content);
  };

  const handleChange = (event: React.ChangeEvent<HTMLTextAreaElement>) => {
    const newValue = event.target.value;
    setSource(newValue);

    if (isValidXYZ(newValue)) {
      setIsValidXYZFormat(true);
    } else {
      setIsValidXYZFormat(false);
    }
  };

  return (
    <Card className={styles.SourceEditor}>
      <VStack>
        <h4>Source editor</h4>
        <ControlGroup>
          <InputGroup readOnly={true} value={sourceName} fill={true} />
          <Tag
            style={{ width: '30%' }}
            intent={isValidXYZFormat ? 'success' : 'warning'}
          >
            {isValidXYZFormat ? 'Correct XYZ pattern' : 'Wrong XYZ pattern'}
          </Tag>

        </ControlGroup>

        <ControlGroup className={styles.ControlGroup}>
          <TextArea
            id='text-editor'
            className={styles.TextEditor}
            value={source}
            onChange={handleChange}
          />
        </ControlGroup>
      </VStack>
    </Card>
  )
}

export default SourceEditor