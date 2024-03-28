import { useLocation } from "react-router-dom";
import PublicationViewPage from "../components/publication-view-page/publication-view-page";
import { publicationService } from "../services/publication/publication.service";
import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import { IPublicationPreview } from "../models/publication-view.model";
import { IVUSSummary } from "../models/vus-summary.model";

const PublicationViewPageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [publications, setPublications] =
    useState<IPublicationPreview[]>(undefined);
  const [variant, setVariant] = useState<IVUSSummary>(undefined);

  const loc = useLocation();
  const variantId = loc.pathname.split("/publication-view/")[1];

  useEffect(() => {
    if (isLoading) {
      publicationService
        ?.getPublicationsByVariantId({ variantId: variantId })
        .then((res) => {
          if (res.isSuccess) {
            setPublications(convertPubDates(res.publications));
            setVariant(res.variant);
            setIsLoading(false);
          } else {
            //TODO: handle error
          }
        });
    }
  }, [isLoading, variantId]);

  if (isLoading) {
    return <Loader />;
  } else {
    return (
      <PublicationViewPage
        variantId={variantId}
        publications={publications}
        variant={variant}
      />
    );
  }

  function convertPubDates(
    publications: IPublicationPreview[]
  ): IPublicationPreview[] {
    let updatedPublications: IPublicationPreview[] = [];

    if (publications) {
      //format publication date
      publications.forEach((publication) => {
        let pub = publication;
        pub.date = new Date(publication.date);
        updatedPublications = updatedPublications.concat(pub);
      });
    }

    return updatedPublications;
  }
};

export default PublicationViewPageWrapper;
