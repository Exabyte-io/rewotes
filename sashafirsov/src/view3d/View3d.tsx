import * as THREE from 'three';
import React, { useCallback, useEffect, useRef, useState } from 'react';
import { Canvas, useFrame, ThreeElements } from '@react-three/fiber';
import { useLocalStorage } from '@uidotdev/usehooks';
import { MeshProps } from '@react-three/fiber';
import { CameraControls } from '@react-three/drei';

import styles from './View3d.module.css';

import Xyz, { XyzSlide } from '../xyz/Xyz';
import { Symbol2Element } from '../Elements';

const DEG45 = Math.PI / 4;


function Box(props: ThreeElements['mesh']) {
    const ref = useRef<THREE.Mesh>(null!);
    const [hovered, hover] = useState(false);
    const [clicked, click] = useState(false);
    return (
        <mesh
            {...props}
            ref={ref}
            scale={clicked ? 1.5 : 1}
            // position={new Vector3(0,0,0)}
            onClick={(event) => click(!clicked)}
            onPointerOver={(event) => hover(true)}
            onPointerOut={(event) => hover(false)}>
            <boxGeometry args={[1, 1, 1, 1, 1, 1]} />
            <meshStandardMaterial color={hovered ? 'hotpink' : 'orange'} wireframe />
        </mesh>
    );
}

function showTooltip(x: number, y: number, text: string) {
    const t = document.getElementById('tooltip')!;
    t.style.left = (x - 380) + 'px';
    t.style.top = (y - 32) + 'px';
    t.style.display = 'block';
    t.style.opacity = '1';
    t.innerText = text;
}

function hideTooltip() {
    const t = document.getElementById('tooltip')!;
    t.style.display = 'none';
}

interface SphereProps extends MeshProps {
    color: string;
    tooltipText: string;
    sourceLine: number;
}

function Sphere(props: SphereProps) {
    const ref = useRef<THREE.Mesh>(null!);
    const [hovered, hover] = useState(false);
    const [clicked, click] = useState(false);
    const [selection] = useLocalStorage<number[]>('EditorSelection');
    const isSelected = props.sourceLine >= selection[0] && props.sourceLine <= selection[1];
    const color = hovered ? 'hotpink' : isSelected ? 'violet' : props.color;
    return (
        <mesh
            {...props}
            ref={ref}
            scale={isSelected ? 2 : clicked ? 1.5 : 1}
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
            <meshStandardMaterial color={color} wireframe />
        </mesh>
    );
}

export default function View3d() {

    const cameraControlRef = useRef<CameraControls | null>(null);

    const [drawing] = useLocalStorage<Xyz>('xyzdrawing');

    const [selection] = useLocalStorage<number[]>('EditorSelection');
    let slide: XyzSlide = new XyzSlide();
    if (drawing && drawing.slides?.length) {
        slide = drawing.slides[0];
        for (const s of drawing.slides)
            if (s.elements.find(e => e.sourceLine >= selection[0] && e.sourceLine <= selection[1])) {
                slide = s;
                break;
            }
    }

    const zoom = 10;
    const keyUpHandler = useCallback((ev: KeyboardEvent) => {

        const key2handler: { [key: string]: () => void } = {
            0: () => cameraControlRef.current?.reset(),
            '-': () => cameraControlRef.current?.zoom(-1, true),
            '_': () => cameraControlRef.current?.zoom(-1, true),
            '=': () => cameraControlRef.current?.zoom(1, true),
            '+': () => cameraControlRef.current?.zoom(1, true),
            'ArrowRight': () => cameraControlRef.current?.rotate(DEG45, 0, true),
            'ArrowLeft': () => cameraControlRef.current?.rotate(-1 * DEG45, 0, true),
            'ArrowUp': () => cameraControlRef.current?.rotate(0, DEG45, true),
            'ArrowDown': () => cameraControlRef.current?.rotate(0, -1 * DEG45, true)
        };
        if (key2handler[ev.key]) {
            ev.preventDefault();
            key2handler[ev.key]();
        }
    }, []);

    useEffect(() => {
        document.addEventListener('keyup', keyUpHandler);
        return () => {
            document.removeEventListener('keyup', keyUpHandler);
        };
    }, [keyUpHandler]);
    return (
        <div className={styles.content}>
            <div id="tooltip" className={styles.tooltip}></div>
            <Canvas camera={{ fov: 75, near: 0.1, far: 1000, position: [zoom, zoom, zoom] }}>
                <CameraControls ref={cameraControlRef} />

                <hemisphereLight position={[0, 1, 0]} args={[0xffffff, 0x888888, 3]} />
                <pointLight position={[10, 10, 10]} />
                <Box position={[-1.2, 0, 0]} />

                {slide.elements.map((e, i) => (
                    <Sphere position={[e.x, e.y, e.z]} key={i}
                            color={Symbol2Element[e.element]?.color ?? 'silver'}
                            sourceLine={e.sourceLine}
                            tooltipText={JSON.stringify(e)}
                    />
                ))}
            </Canvas>
        </div>
    );
}