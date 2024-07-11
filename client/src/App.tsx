import React, { useContext, useEffect, useState } from "react";
import styles from "./App.module.scss";
import {
  Link,
  Route,
  Routes,
  useLocation,
  useNavigate,
} from "react-router-dom";
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
import VusUploadPageWrapper from "./wrappers/vus-upload-page-wrapper";
import ReviewPageWrapper from "./wrappers/review-page-wrapper";
import ReviewHistoryPageWrapper from "./wrappers/review-history-page-wrapper";
import ErrorPage from "./components/error-page/error-page";
import HomePage from "./components/home-page/home-page";
import { AppContext } from "./app-context";
import Banner from "./atoms/banner/banner";
import { IStatus } from "./services/vus/vus.dto";
import Modal from "./atoms/modal/modal";
import VusTable from "./components/view-all-vus-page/vus-table/vus-table";
import Icon from "./atoms/icons/icon";
import HomePageWrapper from "./wrappers/homepage-wrapper";

type AppProps = {};
const App: React.FunctionComponent<AppProps> = () => {
  const navigate = useNavigate();
  let location = useLocation();

  const [isFileUploadModalVisible, setIsFileUploadModalVisible] =
    useState(false);
  const [fileUploadTaskOnDisplay, setFileUploadTaskOnDisplay] =
    useState<IStatus>(undefined);

  const {
    isUserLoggedIn,
    setIsUserLoggedIn,
    completedTasks,
    setCompletedTasks,
  } = useContext(AppContext);

  useEffect(() => {
    authService.isUserLoggedIn().then((res) => {
      setIsUserLoggedIn(res.isUserLoggedIn);

      if (!res.isUserLoggedIn && location.pathname !== "/register") {
        navigate("/login");
      }
    });
  }, [location.pathname]);

  return (
    //TODO: lazy loading
    //routing: https://hygraph.com/blog/routing-in-react

    <div className={styles.container}>
      {completedTasks.length > 0 && (
        <div className={styles["completed-task-ids"]}>
          <div className={styles.banners}>
            {completedTasks.map((t: IStatus) => (
              <div className={styles.banner}>
                <Banner
                  isClosable={true}
                  isGreen={t.isSuccess}
                  onCloseCallback={() =>
                    setCompletedTasks(
                      completedTasks.filter((task) => task.taskId !== t.taskId)
                    )
                  }
                >
                  {t.isSuccess ? (
                    <p>
                      File <b>{t.filename}</b> has been uploaded
                      successfully.&nbsp;
                      <span
                        style={{
                          cursor: "pointer",
                          textDecoration: "underline",
                        }}
                        onClick={() => {
                          setIsFileUploadModalVisible(true);
                          setFileUploadTaskOnDisplay(t);
                        }}
                      >
                        Click here
                      </span>
                      &nbsp; to view upload information.
                    </p>
                  ) : (
                    <p>File {t.filename} failed to upload.</p>
                  )}
                </Banner>
              </div>
            ))}
          </div>
        </div>
      )}

      <Header isUserLoggedIn={isUserLoggedIn} authService={authService} />

      <Routes>
        <Route
          path="/file-upload"
          element={<VusFileUploadPage vusService={vusService} />}
        />
        <Route path="/vus-upload" element={<VusUploadPageWrapper />} />
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
        <Route path="/review/*" element={<ReviewPageWrapper />} />
        <Route
          path="/review-history/*"
          element={<ReviewHistoryPageWrapper />}
        />
        <Route path="/" element={<HomePageWrapper />} />
        <Route path="*" element={<ErrorPage />} />
      </Routes>

      {isFileUploadModalVisible && fileUploadTaskOnDisplay && (
        <Modal
          modalContainerStyle={styles["file-upload-modal"]}
          title={`File Upload Summary`}
          isClosable={true}
          onCloseIconClickCallback={() => setIsFileUploadModalVisible(false)}
        >
          <div className={styles["file-upload-modal-content"]}>
            <p>
              Below you can find a summary of the VUS found in the successfully
              uploaded file&nbsp;
              <b style={{ color: "#008080" }}>
                {fileUploadTaskOnDisplay.filename}
              </b>
              .
            </p>
            {fileUploadTaskOnDisplay.noHpoTermPhenotypes?.length > 0 && (
              <div className={styles["no-hpo-terms-phenotypes-container"]}>
                <div
                  className={styles["no-hpo-terms-phenotypes-title-container"]}
                >
                  <p className={styles["no-hpo-terms-phenotypes-title"]}>
                    <Icon name="warning" fill="#008080" /> Phenotypes warning!
                  </p>
                  <p>
                    No exact match was found with an HPO term for the below
                    phenotypes. Kindly access the sample pages to add existing
                    HPO terms that can replace the inputted phenotypes.
                  </p>
                </div>
                <div className={styles["no-hpo-terms-phenotypes"]}>
                  {fileUploadTaskOnDisplay.noHpoTermPhenotypes.map(
                    (noHpoTermPhenotype) => (
                      <div className={styles["no-hpo-terms-phenotype"]}>
                        <p>
                          <span className={styles.phenotype}>
                            {noHpoTermPhenotype.phenotype}
                          </span>
                          &nbsp;was observed in the following samples:
                        </p>
                        <div
                          className={styles["no-hpo-terms-phenotypes-samples"]}
                        >
                          {noHpoTermPhenotype.samples.map((s) => (
                            // open in new window the sample page
                            <Link
                              to={`/sample/${s}`}
                              target="_blank"
                              rel="noopener noreferrer"
                            >
                              {s}
                            </Link>
                          ))}
                        </div>
                      </div>
                    )
                  )}
                </div>
              </div>
            )}
            Scroll to the right on the table to see which variants' RSIDs need
            to be reviewed. Variants which had been uploaded in the past do not
            get overwritten by the newly uploaded data. For existing variants,
            new samples are added to the variant together with the respective
            HGVS.
            {fileUploadTaskOnDisplay.existingVariantIds.length > 0 && (
              <p>
                <p>
                  These are the variant Ids of the variants which had already
                  been uploaded:
                </p>
                {fileUploadTaskOnDisplay.existingVariantIds.map((id) => (
                  <a
                    href={`/vus/${id}`}
                    target="_blank"
                    rel="noopener noreferrer"
                  >
                    {id},&nbsp;
                  </a>
                ))}
              </p>
            )}
            {fileUploadTaskOnDisplay.vusList && (
              <div className={styles["file-content"]}>
                <div className={styles["file-summary"]}>
                  <div className={styles.summary}>
                    <p>RSIDs</p>
                    <p>
                      {
                        fileUploadTaskOnDisplay.vusList.filter(
                          (vus) =>
                            vus.rsid.length > 0 &&
                            vus.rsid !== "NORSID" &&
                            !vus.rsidReviewRequired
                        ).length
                      }
                    </p>
                  </div>
                  <div className={`${styles.summary} ${styles["num-vus"]}`}>
                    <p>VUS</p> <p>{fileUploadTaskOnDisplay.vusList.length}</p>
                  </div>
                  <div className={styles.summary}>
                    <p>ClinVar</p>
                    <p>
                      {
                        fileUploadTaskOnDisplay.vusList.filter(
                          (vus) => vus.isFoundInClinvar
                        ).length
                      }
                    </p>
                  </div>
                </div>

                <VusTable
                  vusList={fileUploadTaskOnDisplay.vusList}
                  onVariantClickCallback={() =>
                    setIsFileUploadModalVisible(false)
                  }
                />
              </div>
            )}
          </div>
        </Modal>
      )}
    </div>
  );
};

export default App;
