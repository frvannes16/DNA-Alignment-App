import Console from './components/Console.tsx';
import * as React from 'react';
import * as ReactDOM from 'react-dom';

const wrapper = document.getElementById('app');
wrapper ? ReactDOM.render(<Console />, wrapper) : null;