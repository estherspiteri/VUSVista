const dotenv = require('dotenv');
const environment = process.env.NODE_ENV || 'development';
const envFile = `.env.${environment}`;

console.log(`Loading environment variables from ${envFile}`);
dotenv.config({ path: envFile });