# Pull official base image
FROM node:18-alpine AS build

# Set working directory inside the container
WORKDIR /client

# Add `/client/node_modules/.bin` to $PATH
ENV PATH /client/node_modules/.bin:$PATH

# Copy the package.json and package-lock.json to the container
COPY package.json package-lock.json ./

# Copy the rest of the application files into the container
COPY . .

# Install dependencies
RUN npm ci
RUN npm run build

# Use an official Nginx runtime as a parent image
FROM nginx:stable-alpine

# Copy the ngnix.conf to the container
COPY --from=build /client/nginx.conf /etc/nginx/conf.d/default.conf

# Copy the React app build files to the container
COPY --from=build /client/build /usr/share/nginx/html

# Expose port 80 for Nginx
EXPOSE 80

# Start Nginx when the container starts
CMD ["nginx", "-g", "daemon off;"]

