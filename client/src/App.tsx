import React, { useEffect, useState } from "react";
import styles from "./App.module.scss";
import { Route, Routes, useLocation, useNavigate } from "react-router-dom";
import { vusService } from "./services/vus/vus.service";
import Header from "./components/header/header";
import ViewVusPage from "./components/view-all-vus-page/view-all-vus-page";
import VusFileUploadPage from "./components/vus-file-upload-page/vus-file-upload-page";
import ViewAllSamples from "./components/view-samples-page/view-samples-page";
import { samplesService } from "./services/sample/sample.service";
import PublicationViewPageWrapper from "./wrappers/publication-view-page-wrapper";
import LoginPage from "./components/login-page/login-page";
import { authService } from "./services/auth/auth.service";
import RegisterPage from "./components/register-page/register-page";
import ProfilePageWrapper from "./wrappers/profile-page-wrapper";
import SamplePageWrapper from "./wrappers/sample-page-wrapper";
import VusPageWrapper from "./wrappers/vus-page-wrapper";
import PublicationPhenotypeViewPageWrapper from "./wrappers/publication-phenotype-view-page-wrapper";
import VusUploadPage from "./components/vus-upload-page/vus-upload-page";

type AppProps = {};
//TODO: add session cookie Id
const App: React.FunctionComponent<AppProps> = () => {
  const [isUserLoggedIn, setIsUserLoggedIn] = useState(false);
  const navigate = useNavigate();
  let location = useLocation();

  //TODO: is this correct location?
  useEffect(() => {
    authService.isUserLoggedIn().then((res) => {
      setIsUserLoggedIn(res.isUserLoggedIn);

      if (!res.isUserLoggedIn && location.pathname !== "/register") {
        navigate("/login");
      } else if (location.pathname !== "/login") {
        //TODO: redirect to profile page
      }
    });
  }, [location.pathname]);

  return (
    //TODO: lazy loading
    //routing: https://hygraph.com/blog/routing-in-react
    <div className={styles.container}>
      <Header isUserLoggedIn={isUserLoggedIn} authService={authService} />
      <Routes>
        <Route
          path="/file-upload"
          element={<VusFileUploadPage vusService={vusService} />}
        />
        <Route
          path="/vus-upload"
          element={
            <VusUploadPage
              sampleService={samplesService}
              vusService={vusService}
            />
          }
        />
        <Route
          path="/view-samples"
          element={<ViewAllSamples sampleService={samplesService} />}
        />
        <Route path="/sample/*" element={<SamplePageWrapper />} />
        <Route
          path="/view-vus"
          element={<ViewVusPage vusService={vusService} />}
        />
        <Route path="/vus/*" element={<VusPageWrapper />} />
        <Route
          path="/publication-view/*"
          element={<PublicationViewPageWrapper />}
        />
        <Route
          path="/publication-phenotype-view/*"
          element={<PublicationPhenotypeViewPageWrapper />}
        />
        <Route
          path="/login"
          element={<LoginPage authService={authService} />}
        />
        <Route
          path="/register"
          element={<RegisterPage authService={authService} />}
        />
        <Route path="/profile" element={<ProfilePageWrapper />} />
        {/*TODO: handle no route match*/}
        {/* <Route path="*" element={<NoMatch />} /> */}
      </Routes>
    </div>
  );
};

export default App;
