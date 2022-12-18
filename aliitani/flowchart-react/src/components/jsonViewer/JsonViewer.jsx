import { JsonDataState } from '../../state/AppState';
import './jsonViewer.css';
import { useSnapshot } from 'valtio';

const JsonViewer = () => {
	const jsonState = useSnapshot(JsonDataState);
	return (
		<div id='json-viewer' className='overflow-scroll overflow-x-hidden'>
			<p>JSON Data:</p>

			<pre>{JSON.stringify(jsonState, undefined, 2)}</pre>
		</div>
	);
};

export default JsonViewer;
