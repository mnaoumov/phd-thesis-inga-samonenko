using System.Web.Mvc;

namespace PhdThesisIngaSamonenko.Controllers
{
    public class KMDController : Controller
    {
        public ActionResult Index()
        {
            return File("KMD.R", "text/plain");
        }

    }
}
