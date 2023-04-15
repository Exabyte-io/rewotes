import { Button, Card, ControlGroup, Tag } from '@blueprintjs/core'
import React, { useEffect, useRef, useState } from 'react'
import styles from './StructureViewer.module.scss';

import { Object3D, Vector3 } from 'three';
import { Canvas, extend, useFrame, useThree } from '@react-three/fiber';
import { Sphere,  } from '@react-three/drei';

import ViewHeading from '../../components/view_heading/ViewHeading';
import { useSettings } from '../../context/SettingsContext';
import VStack from '../../components/utils/VStack';
import { useAtomsContext, Atom } from '../../context/AtomsContext';
import { useSourceContext } from '../../context/SourceContext';

import { DragControls } from 'three/examples/jsm/controls/DragControls';
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";
extend({ OrbitControls });

function useOrbitControls(
  ref: React.MutableRefObject<OrbitControls | null>,
  camera: THREE.Camera,
  domElement: HTMLElement,
) {
  useEffect(() => {
    const controls = new OrbitControls(camera, domElement);
    ref.current = controls;
    return () => {
      controls.dispose();
    };
  }, [ref, camera, domElement]);
}
interface ControlsProps {
  shiftPressed: boolean;
  editEnabled: boolean | undefined;
}

function Controls(props: ControlsProps) {
  const { camera, gl } = useThree();
  const ref = useRef<OrbitControls | null>(null);
  useOrbitControls(ref, camera, gl.domElement);

  useEffect(() => {
    if (ref.current) {
      ref.current.enabled = !(props.editEnabled && props.shiftPressed);
    }
  }, [props.editEnabled, props.shiftPressed]);

  return null;
}

interface DraggableSphereProps {
  atom: Atom;
  materialColor: string | undefined;
  onAtomMove: (atom: Atom, newPosition: Vector3) => void;
  editEnabled?: boolean;
  shiftPressed?: boolean;
};

const DraggableSphere: React.FC<DraggableSphereProps> = (props) => {
  const { atom, materialColor, onAtomMove, editEnabled, shiftPressed } = props;
  const mesh = useRef<THREE.Mesh>(null);
  const { camera, gl } = useThree();

  useEffect(() => {
    if (!editEnabled || !mesh.current || !shiftPressed) return;

    const controls = new DragControls([mesh.current], camera, gl.domElement);

    controls.addEventListener('drag', (e) => {
      console.log('drag');
      if (!shiftPressed) {
        e.preventDefault();
        e.stopPropagation();
        return false;
      }
      
      const x = e.object.position.x;
      const y = e.object.position.y;
      const z = e.object.position.z;
      atom.position.set(x, y, z);
    });

    controls.addEventListener('dragend', (e) => {
      onAtomMove(atom, atom.position);
    });

    return () => controls.dispose();
  }, [camera, gl, atom, onAtomMove, editEnabled, shiftPressed]);

  return (
    <Sphere args={[0.1]} position={atom.position.toArray()} ref={mesh}>
      <meshStandardMaterial attach="material" color={materialColor} />
    </Sphere>
  );
};

const StructureViewer: React.FC = () => {
  const settings = useSettings();
  const theme = settings.settings.theme;
  const editEnabled = settings.settings.editingIn3D;
  const atomsDisplayData = settings.settings.atomsDisplayData || [];

  const viewBoxSize: number = 20;
  const {atoms} = useAtomsContext();
  const {source, setSource} = useSourceContext();
  const [shiftPressed, setShiftPressed] = useState(false);

  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Shift') setShiftPressed(true);
    };
  
    const handleKeyUp = (e: KeyboardEvent) => {
      if (e.key === 'Shift') setShiftPressed(false);
    };
  
    window.addEventListener('keydown', handleKeyDown);
    window.addEventListener('keyup', handleKeyUp);
  
    return () => {
      window.removeEventListener('keydown', handleKeyDown);
      window.removeEventListener('keyup', handleKeyUp);
    };
  }, []);

  
  Object3D.DEFAULT_UP.set(0, 0, 1); // By scientific convention Z is up: (0,0,1), by default ThreeJS uses (0,1,0) as up
  // TODO: set camera default position facing (1,1,0)


  const handleAtomMove = (atom: Atom, newPosition: Vector3) => {
    const atomIndex = atoms.findIndex((a) => a.id === atom.id);
    if (atomIndex === -1) return;
  
    const sourceLines = source.split('\n');
    const updatedLine = `${atoms[atomIndex].element} ${newPosition.x.toFixed(3)} ${newPosition.y.toFixed(3)} ${newPosition.z.toFixed(3)}`;
    sourceLines[atomIndex + 2] = updatedLine;
    const newSourceData = sourceLines.join('\n');
    setSource(newSourceData);
    
  };

  function toggleEditMode() {
    settings.updateSettings({editingIn3D: !editEnabled})
  }

  return (
    <Card className={styles.StructureViewer}>
      <VStack>
      <ViewHeading>
        <h4>Structure Viewer</h4>
      </ViewHeading>
      <ControlGroup>
        <Button onClick={toggleEditMode} active={editEnabled}>
          Edit 3D
          </Button>
          {editEnabled && <Tag>Press Shift to move atoms</Tag>}
      </ControlGroup>
    <Canvas className={styles.Canvas + ' ' + theme}>
      <Controls  shiftPressed={shiftPressed} editEnabled={editEnabled} />
      <ambientLight />
      <pointLight position={[viewBoxSize, viewBoxSize, viewBoxSize]} />
      <gridHelper args={[viewBoxSize, viewBoxSize, '#CCCCCC', '#AAAAAA']} rotation-x={Math.PI / 2} />
      <axesHelper args={[viewBoxSize/2]} />
      {atoms.map((atom, index) => {
            const colorObj = atomsDisplayData.find((color) => color.element === atom.element);
            const defaultColor = settings.settings.defaultAtomColor;
            const materialColor = colorObj ? colorObj.color : defaultColor;
            return <DraggableSphere 
              atom={atom} 
              materialColor={materialColor} 
              key={atom.id} 
              onAtomMove={handleAtomMove}
              editEnabled={editEnabled}
              shiftPressed={shiftPressed}
              />;
          })}
      </Canvas>
      </VStack>
    </Card>
  );
};

export default StructureViewer