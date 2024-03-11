import React from "react";
import styles from "./profile-page.module.scss";
import { IProfile } from "../../models/profile.model";

type ProfilePageProps = { profile: IProfile };

//TODO: update profile & update password
const ProfilePage: React.FunctionComponent<ProfilePageProps> = (
  props: ProfilePageProps
) => {
  return (
    <div className={styles["profile-page-container"]}>
      <div className={styles.title}>Profile</div>

      <div className={styles["profile-content"]}>
        <div className={styles["field-container"]}>
          {/** Name */}
          <div className={styles.field}>
            <span className={styles["field-name"]}>Name:</span>
            <span>{props.profile.name}</span>
          </div>

          {/** Surname */}
          <div className={styles.field}>
            <span className={styles["field-name"]}>Surname:</span>
            <span>{props.profile.surname}</span>
          </div>

          {/** Email - TODO: add validation */}
          <div className={styles.field}>
            <span className={styles["field-name"]}>Email:</span>
            <span>{props.profile.email}</span>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ProfilePage;
