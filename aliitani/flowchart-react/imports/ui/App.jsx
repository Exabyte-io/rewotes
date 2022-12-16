import React from 'react';
import { FlowChart } from './FlowChart';
import { JsonViewer } from './JsonViewer';

export const App = () => (
	<div className='max-w-3xl min-h-screen mx-auto sm:pt-10'>
		<FlowChart />
		<JsonViewer />
	</div>
);
