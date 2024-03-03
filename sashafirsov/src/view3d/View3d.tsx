import * as THREE from 'three';
import React, { useRef, useState } from 'react';
import { Canvas, useFrame, ThreeElements } from '@react-three/fiber';
import { useLocalStorage } from '@uidotdev/usehooks';
import { MeshProps } from '@react-three/fiber';

import styles from './View3d.module.css';

import Xyz from '../xyz/Xyz';
import { Symbol2Element } from '../Elements';

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

interface SphereProps extends MeshProps {
    color: string;
    tooltipText: string;
}

function showTooltip(x: number, y: number, text: string) {
    const t = document.getElementById('tooltip')!;
    console.log({ x, y, text });
    t.style.left = x + 'px';
    t.style.top = (y - 32) + 'px';
    t.style.display = 'block';
    t.style.opacity = '1';
    t.innerText = text;
}

function hideTooltip() {
    const t = document.getElementById('tooltip')!;
    t.style.display = 'none';
}

function Sphere(props: SphereProps) {
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
            onPointerOver={(event) => {
                showTooltip(event.x, event.y, props.tooltipText);
                hover(true);
            }}
            onPointerOut={(event) => {
                hideTooltip();
                hover(false);
            }}
            userData={{ tooltipText: props.tooltipText }}
        >
            <sphereGeometry args={[0.1, 32, 32]} />
            <meshStandardMaterial color={hovered ? 'hotpink' : props.color} />
        </mesh>
    );
}

export default function View3d() {
    const [drawing, saveDrawing] = useLocalStorage<Xyz>('xyzdrawing');


    return (
        <div className={styles.content}>
            <div id="tooltip" className={styles.tooltip}></div>
            <Canvas camera={{ fov: 75, near: 0.1, far: 1000, position: [5, 5, 5] }}>
                <hemisphereLight position={[0, 1, 0]} args={[0xffffff, 0x888888, 3]} />
                <pointLight position={[10, 10, 10]} />
                <Box position={[-1.2, 0, 0]} />

                {drawing && drawing.slides[0]?.elements.map((e, i) => (
                    <Sphere position={[e.x, e.y, e.z]} key={i}
                            color={Symbol2Element[e.element]?.color ?? 'silver'}
                            tooltipText={JSON.stringify(e)}
                    />
                ))}
            </Canvas>
        </div>
    );
}