import Console from './components/Console.jsx';
import React from 'react';
import ReactDOM from 'react-dom';

const wrapper = document.getElementById('app');
wrapper ? ReactDOM.render(<Console />, wrapper) : null;