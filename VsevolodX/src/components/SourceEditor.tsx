import { Button, Card, ControlGroup, Dialog, FileInput, TextArea, InputGroup, Tag } from '@blueprintjs/core'
import React, { useState } from 'react'
import styles from './SourceEditor.module.scss';
import { useSourceContext } from './SourceContext';
import { useSettings } from './SettingsContext';

function SourceEditor() {
  const { source, setSource, importSource, sourceName, setSourceName, isValidXYZFormat, setIsValidXYZFormat } = useSourceContext();
  const [isDialogOpen, setIsDialogOpen] = useState(false);

  const theme = useSettings().settings.theme;

  const handleImport = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      importSource(file);
      setSourceName(file.name);
    }
  };

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
      <div className='VStack'>
        
      <h4>Source editor</h4>
      <ControlGroup>
      <Button 
        icon='insert'
        onClick={() => setIsDialogOpen(true)}>Import file</Button>
      <Dialog className={`bp4-${theme}`} isOpen={isDialogOpen} onClose={() => setIsDialogOpen(false)} title="Import file">
        <FileInput text="Choose file..." onInputChange={handleImport} />
      </Dialog>

      <Button icon='import'>Export file</Button>
      </ControlGroup>
      <ControlGroup>
      <InputGroup readOnly={true} value={sourceName} />
      
      <Tag
        title={isValidXYZFormat? 'Correct' : 'Wrong XYZ pattern'}
        intent={isValidXYZFormat? 'success' : 'warning'}
      >
        {isValidXYZFormat? 'Correct' : 'Wrong XYZ pattern'}
      </Tag>

      </ControlGroup>

      <ControlGroup className={styles.ControlGroup}>
          <TextArea  
            className={styles.TextEditor}
            value={source} 
            onChange={handleChange}
            />
      </ControlGroup>
            </div>
    </Card>
  )
}

export default SourceEditor