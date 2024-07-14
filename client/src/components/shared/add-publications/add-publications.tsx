import React, { useState } from "react";
import styles from "./add-publications.module.scss";
import Icon from "../../../atoms/icons/icon";
import Text from "../../../atoms/text/text";

type AddPublications = {
  isAddingPublications?: boolean;
  onPublicationsUpdateCallback?: (publicationUrls: string[]) => void;
};

const AddPublications: React.FunctionComponent<AddPublications> = (
  props: AddPublications
) => {
  const [url, setUrl] = useState("");
  const [addedUrls, setAddedUrls] = useState<string[]>([]);

  return (
    <div className={styles["add-publications-container"]}>
      <div className={styles["new-input"]}>
        <div className={styles.text}>
          <Text
            disabled={props.isAddingPublications}
            placeholder="Insert URL ..."
            value={url}
            autoFocus={true}
            onChange={(e) => setUrl(e.currentTarget.value)}
          />
        </div>
        <div className={styles["icon-wrapper"]}>
          <Icon
            name="add"
            className={styles.add}
            onClick={() => {
              if (!props.isAddingPublications && url.trim().length > 0) {
                if (!addedUrls.find((u) => u === url)) {
                  const updatedUrls = addedUrls.concat(url);
                  setAddedUrls(updatedUrls);

                  props.onPublicationsUpdateCallback &&
                    props.onPublicationsUpdateCallback(updatedUrls);
                }
                setUrl("");
              }
            }}
          />
        </div>
      </div>
      {addedUrls.map((u) => (
        <div className={styles["added-url"]}>
          <a href={u}>{u}</a>
          <Icon
            name="close"
            stroke="#008080"
            onClick={() => {
              const updatedUrls = addedUrls.filter((url) => url !== u);
              setAddedUrls(updatedUrls);

              props.onPublicationsUpdateCallback &&
                props.onPublicationsUpdateCallback(updatedUrls);
            }}
          />
        </div>
      ))}
    </div>
  );
};

export default AddPublications;
