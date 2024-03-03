import * as THREE from 'three';
import React, { useRef, useState } from 'react';
import { Canvas, useFrame, ThreeElements } from '@react-three/fiber';
import { useLocalStorage } from '@uidotdev/usehooks';

import styles from './View3d.module.css';

import Xyz, { ElementXyz, XyzSlide } from '../xyz/Xyz';

function Box(props: ThreeElements['mesh']) {
    const ref = useRef<THREE.Mesh>(null!);
    const [hovered, hover] = useState(false);
    const [clicked, click] = useState(false);
    useFrame((state, delta) => (ref.current.rotation.x += delta));
    return (
        <mesh
            {...props}
            ref={ref}
            scale={clicked ? 1.5 : 1}
            onClick={(event) => click(!clicked)}
            onPointerOver={(event) => hover(true)}
            onPointerOut={(event) => hover(false)}>
            <boxGeometry args={[1, 1, 1]} />
            <meshStandardMaterial color={hovered ? 'hotpink' : 'orange'} />
        </mesh>
    );
}


function Sphere(props: ThreeElements['mesh']) {
    const ref = useRef<THREE.Mesh>(null!);
    const [hovered, hover] = useState(false);
    const [clicked, click] = useState(false);
    useFrame((state, delta) => (ref.current.rotation.x += delta));
    return (
        <mesh
            {...props}
            ref={ref}
            scale={clicked ? 1.5 : 1}
            onClick={(event) => click(!clicked)}
            onPointerOver={(event) => hover(true)}
            onPointerOut={(event) => hover(false)}>
            <sphereGeometry args={[0.1, 32, 32]} />
            <meshStandardMaterial color={hovered ? 'hotpink' : 'orange'} />
        </mesh>
    );
}

export default function View3d() {
    const [drawing, saveDrawing] = useLocalStorage<Xyz>('xyzdrawing');


    return (
        <div className={styles.content}>
            <Canvas camera={{ fov: 75, near: 0.1, far: 1000, position: [5, 5, 5] }}>
                <hemisphereLight position={[0, 1, 0 ]} args={[0xffffff, 0x888888, 3 ]}/>
                <pointLight position={[10, 10, 10]} />
                <Box position={[-1.2, 0, 0]} />
                <Sphere position={[1.2, 0, 0]} />

                {drawing&& drawing.slides[0].elements.map((e,i)=>(
                    <Sphere position={[e.x, e.y, e.z]} key={i}/>
                ))}
            </Canvas>
        </div>
    );
}