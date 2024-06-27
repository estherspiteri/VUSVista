const customFetch = async (url, options) => {
  try {
    const response = await fetch(url, options);
    console.log("response", response);
    if (!response.ok) {
      // Handle HTTP errors by navigating to the error page
      window.location.href = "/error";
    }
    
    // Attempt to parse the JSON response
    try {
      const resJson = await response.json();
      return resJson;
    } catch (jsonError) {
      console.error("JSON Parsing Error:", jsonError);
      window.location.href = "/error";
    }
  } catch (error) {
    // Handle fetch errors (e.g., network issues)
    console.error("Fetch Error:", error);
    window.location.href = "/error";
  }
};

export default customFetch;
